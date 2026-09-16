"""The parent process's read-only view of the reference data.

The app needs two things from the bundle before it schedules anything: whether a
transcript exists, and the release and ClinVar versions that go into the cache key.
Worker processes open their own sources; nothing here is shared across processes.

Connections are per-thread. The ASGI server answers requests on a thread pool, and a
single sqlite3 connection used concurrently interleaves its cursors: one thread steps
another's statement, `fetchone` comes back empty, and a transcript that is in the
bundle is reported as missing. Read-only connections cost about 0.2 ms to open and
the pool is bounded, so one per thread is the cheap fix rather than a lock that would
serialise every lookup.
"""
import os
import re
import sqlite3
import threading

from bedesign.transcript_source import (TRANSCRIPTS_DB, find_bundle,
                                        find_clinvar_db)

CLINVAR_VERSION = re.compile(r'clinvar-(.+)\.db$')

# Enough to show that a search was too broad without rendering a wall of symbols.
GENE_LIMIT = 25
# The busiest human genes have a few hundred transcripts; MAP2K1 has 32.
TRANSCRIPT_LIMIT = 200
TRANSCRIPT_FIELDS = ('transcript_id', 'display_name', 'biotype', 'seq_region',
                     'strand', 'start', 'end', 'mane_select', 'ensembl_canonical',
                     'cds_length')


class References(object):
    def __init__(self, refdata='refdata', clinvar_db=''):
        self.bundle = find_bundle(refdata)
        self.clinvar_db = find_clinvar_db(clinvar_db or None, refdata)
        self._uri = 'file:%s?mode=ro' % os.path.join(self.bundle, TRANSCRIPTS_DB)
        self._local = threading.local()
        meta = dict(self._db().execute('SELECT key, value FROM meta').fetchall())
        self.release = meta.get('release', 'unknown')
        match = CLINVAR_VERSION.search(os.path.basename(self.clinvar_db))
        self.clinvar_version = match.group(1) if match else 'unknown'

    def _db(self):
        """This thread's connection, opened on first use."""
        db = getattr(self._local, 'db', None)
        if db is None:
            db = sqlite3.connect(self._uri, uri=True)
            self._local.db = db
        return db

    def transcript_exists(self, transcript_id):
        """Parameterized, and only ever called with an ID that already matched
        the transcript pattern."""
        row = self._db().execute(
            'SELECT 1 FROM transcript WHERE transcript_id = ? LIMIT 1',
            (transcript_id,)).fetchone()
        return row is not None

    def search_genes(self, prefix, limit=GENE_LIMIT):
        """Gene symbols starting with `prefix`, for the search page.

        A range predicate, not LIKE. `idx_transcript_gene_name` is a BINARY-collation
        index and SQLite's LIKE is case-insensitive by default, so LIKE cannot use it:
        measured on the 646,577-row bundle, `LIKE 'MAP2K%'` scans the table in 10.7 ms
        while this is index-backed at 0.04 ms. On a search box that is a 250x
        amplification factor kept off the abuse surface, not just a speedup.

        The upper bound is the prefix with the highest code point appended, which is
        what makes '>= prefix AND < bound' mean 'starts with prefix'.
        """
        if not prefix:
            return []
        rows = self._db().execute(
            'SELECT DISTINCT gene_name FROM transcript '
            'WHERE gene_name >= ? AND gene_name < ? '
            'ORDER BY gene_name LIMIT ?',
            (prefix, prefix + '\uffff', limit)).fetchall()
        return [row[0] for row in rows]

    def resolve_gene(self, symbol):
        """The gene a symbol names, and every symbol it could have named.

        An exact hit wins over a prefix list, and a prefix that can only mean one gene
        resolves too, so searching 'MAP2K1' does not stop to ask which MAP2K1. An
        empty gene means the caller should show the matches and let the user choose.
        """
        matches = self.search_genes(symbol)
        if symbol in matches:
            return symbol, matches
        if len(matches) == 1:
            return matches[0], matches
        return '', matches

    def transcripts_for_gene(self, gene_name, limit=TRANSCRIPT_LIMIT):
        """Every transcript of a gene, best choice first.

        MANE Select is the transcript RefSeq and Ensembl agree on, so it leads;
        Ensembl canonical is the fallback for genes that have none, and protein-coding
        sorts above the retained-intron and NMD entries that a base editor screen
        rarely wants.

        The projection is explicit because `cds_sequence` and `protein_sequence` are
        why this database is 1.07 GB -- SELECT * here would read a megabyte to render
        a list.
        """
        rows = self._db().execute(
            'SELECT transcript_id, display_name, biotype, seq_region, strand, '
            '       start, end, mane_select, ensembl_canonical, '
            '       length(cds_sequence) '
            'FROM transcript WHERE gene_name = ? '
            'ORDER BY mane_select DESC, ensembl_canonical DESC, '
            "         biotype = 'protein_coding' DESC, display_name LIMIT ?",
            (gene_name, limit)).fetchall()
        return [dict(zip(TRANSCRIPT_FIELDS, row)) for row in rows]
