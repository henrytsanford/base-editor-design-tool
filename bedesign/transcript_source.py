"""Where reference data comes from.

The design script needs five things from Ensembl per transcript. This module puts
them behind one interface, served either by rest.ensembl.org (`EnsemblRestSource`)
or by a local reference bundle built by tools/build_reference.py (`LocalSource`).

It also serves the other reference set a run needs, the ClinVar variants for a gene
(`ClinVarSource`), out of a database built by tools/build_clinvar.py.
"""
import glob
import json
import os
import random
import re
import sqlite3
import time

import pandas as pd
import requests
from Bio.Seq import reverse_complement

ENSEMBL_SERVER = "https://rest.ensembl.org"
ENSEMBL_ATTEMPTS = 7
ENSEMBL_TIMEOUT = 90
ENSEMBL_MAX_BACKOFF = 60
_session = requests.Session()

JSON_HEADERS = {"Content-Type": "application/json"}
TEXT_HEADERS = {"Content-Type": "text/plain"}

# Everything a bundle must contain to serve a design run.
TRANSCRIPTS_DB = 'transcripts.db'
GENOME_FA = 'genome.fa.bgz'
BUNDLE_FILES = (TRANSCRIPTS_DB, GENOME_FA)

# The ClinVar database written by tools/build_clinvar.py. The column names are the
# internal ones the design script is written against, not ClinVar's own.
# Columns tools/build_reference.py writes, and the ones a bundle must have for
# LocalSource to answer a lookup. Shared so writer, reader and tests agree.
TRANSCRIPT_TABLE = 'transcript'
TRANSCRIPT_COLUMNS = (
    'transcript_id', 'display_name', 'gene_id', 'gene_name', 'biotype',
    'seq_region', 'strand', 'start', 'end', 'exons', 'cds', 'cds_sequence',
    'protein_sequence', 'mane_select', 'ensembl_canonical',
)
REQUIRED_TRANSCRIPT_COLUMNS = ('mane_select', 'ensembl_canonical')

CLINVAR_DB_GLOB = 'clinvar-*.db'
VARIANT_TABLE = 'variant'
VARIANT_COLUMNS = ('#AlleleID', 'RefSeqID', 'Name', 'GeneSymbol',
                   'ClinicalSignificance', 'PhenotypeList', 'ClinVar_SNP_Position',
                   'ReferenceAllele', 'AlternateAllele', 'ReviewStatus')
# One gene is queried once per edit in a run, so a handful of symbols is plenty.
VARIANT_CACHE_SIZE = 4


class EnsemblUnavailable(Exception):
    """The reference service could not be reached after repeated attempts."""


class BundleNotFound(Exception):
    """No usable local reference bundle was found on disk."""


class ClinVarNotFound(Exception):
    """No usable ClinVar database was found on disk."""


class TranscriptNotFound(Exception):
    """The transcript ID is not in the reference data."""

    def __init__(self, tr, detail=''):
        Exception.__init__(self,
            "Transcript '%s' not found in Ensembl%s.\n"
            "Check the ID at https://www.ensembl.org; it should be a bare transcript "
            "ID such as ENST00000294952 (version suffixes like .13 are stripped "
            "automatically)." % (tr, ' (%s)' % detail if detail else ''))


'''
GETs an Ensembl REST endpoint, retrying transient failures. rest.ensembl.org
regularly returns 500/503 under load, in correlated bursts, and a healthy response
can take 20s. Non-retryable responses (including 4xx) are returned so callers can
decide what they mean; a persistent outage raises EnsemblUnavailable.
'''
def ensembl_get(ext, headers, attempts=ENSEMBL_ATTEMPTS, timeout=ENSEMBL_TIMEOUT):
    last = None
    for attempt in range(1, attempts + 1):
        try:
            r = _session.get(ENSEMBL_SERVER + ext, headers=headers, timeout=timeout)
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            last = type(e).__name__
            retry_after = None
        else:
            if r.status_code != 429 and r.status_code < 500:
                return r
            last = 'HTTP ' + str(r.status_code)
            retry_after = r.headers.get('Retry-After')
        if attempt == attempts:
            break
        if retry_after:
            try:
                delay = float(retry_after)
            except ValueError:
                delay = 2 ** attempt
        else:
            delay = 2 ** attempt + random.uniform(0, 1)
        # cap the wait even if Retry-After asks for something absurd
        delay = min(delay, ENSEMBL_MAX_BACKOFF)
        print('Ensembl request failed (%s), retrying in %.0fs (attempt %d of %d)'
              % (last, delay, attempt + 1, attempts))
        time.sleep(delay)
    raise EnsemblUnavailable(
        "Ensembl REST API is unavailable: %s after %d attempts for %s"
        % (last, attempts, ext)
    )


class TranscriptSource(object):
    """The five pieces of reference data a design run needs.

    Return shapes match the Ensembl REST responses the script was written
    against, so implementations are interchangeable at every call site.
    """

    def describe(self):
        """One line naming the reference data, recorded in the run README."""
        raise NotImplementedError

    def lookup(self, tr):
        """Transcript record: /lookup/id/{tr}?expand=1 JSON.

        Uses display_name, assembly_name, strand, seq_region_name, Parent,
        start, end, and Exon[] (each with start/end).
        Raises TranscriptNotFound if the ID is unknown.
        """
        raise NotImplementedError

    def cds_mappings(self, tr, length):
        """CDS blocks in genomic coordinates: the 'mappings' list of /map/cds.
        None when the transcript has no CDS (e.g. non-coding)."""
        raise NotImplementedError

    def genomic_sequence(self, tr, flank=40):
        """Genomic span of the transcript in transcript orientation, with
        `flank` bp added at each end. '' if unavailable."""
        raise NotImplementedError

    def protein_sequence(self, tr):
        """Translated protein sequence. '' if unavailable."""
        raise NotImplementedError

    def cds_sequence(self, tr):
        """Spliced coding sequence. '' if unavailable."""
        raise NotImplementedError


class EnsemblRestSource(TranscriptSource):
    """Live rest.ensembl.org, with the retry policy in `ensembl_get`."""

    def describe(self):
        return 'Ensembl REST API (%s)' % ENSEMBL_SERVER

    def lookup(self, tr):
        r = ensembl_get("/lookup/id/" + tr + "?expand=1", JSON_HEADERS)
        if not r.ok:
            raise TranscriptNotFound(tr, 'HTTP %d' % r.status_code)
        return r.json()

    def cds_mappings(self, tr, length):
        try:
            r = ensembl_get("/map/cds/" + tr + "/1.." + str(length) + "?", JSON_HEADERS)
        except EnsemblUnavailable:
            # /map/cds answers a non-coding transcript with a persistent 500, which
            # by status alone looks like an outage (confirmed live: MALAT1 500s
            # while a coding transcript 200s). If the API is answering for this
            # transcript, the 500 means "no CDS" rather than "down".
            if not self.reachable(tr):
                raise
            return None
        if not r.ok:
            return None
        # get_exons indexes [0], so an empty mapping list must come back as None
        return r.json()['mappings'] or None

    def reachable(self, tr):
        """Is the API serving this transcript right now? One attempt, no retries."""
        try:
            return ensembl_get("/lookup/id/" + tr, JSON_HEADERS, attempts=1).ok
        except EnsemblUnavailable:
            return False

    def genomic_sequence(self, tr, flank=40):
        r = ensembl_get(
            "/sequence/id/" + tr + "?content-type=text/plain;expand_5prime="
            + str(flank) + ";expand_3prime=" + str(flank), TEXT_HEADERS)
        return r.text if r.ok else ''

    def protein_sequence(self, tr):
        r = ensembl_get("/sequence/id/" + tr + "?content-type=text/plain;type=protein",
                        TEXT_HEADERS)
        return r.text if r.ok else ''

    def cds_sequence(self, tr):
        r = ensembl_get("/sequence/id/" + tr + "?content-type=text/plain;type=cds",
                        TEXT_HEADERS)
        return r.text if r.ok else ''


def missing_bundle_files(path):
    return [f for f in BUNDLE_FILES if not os.path.exists(os.path.join(path, f))]


class LocalSource(TranscriptSource):
    """A reference bundle on disk, built by tools/build_reference.py.

    Answers the same five questions as the REST API with no network, pinned to
    one Ensembl release so a rerun gives the same guides.
    """

    def __init__(self, bundle):
        missing = missing_bundle_files(bundle)
        if missing:
            raise BundleNotFound(
                "Reference bundle at '%s' is missing %s.\n"
                "Build one with: python tools/build_reference.py"
                % (bundle, ', '.join(missing)))
        self.bundle = bundle
        self.genome_path = os.path.join(bundle, GENOME_FA)
        self._db = sqlite3.connect(
            'file:%s?mode=ro' % os.path.join(bundle, TRANSCRIPTS_DB), uri=True)
        self._db.row_factory = sqlite3.Row
        have = {c['name'] for c in
                self._db.execute('PRAGMA table_info(transcript)')}
        stale = [c for c in REQUIRED_TRANSCRIPT_COLUMNS if c not in have]
        if stale:
            # A bundle is a build artifact, not user data, so an old one is
            # rebuilt rather than worked around. Tolerating the missing columns
            # would make every gene look like it has no MANE transcript, which
            # reads as a bug in the caller rather than a stale bundle.
            raise BundleNotFound(
                "Reference bundle at '%s' predates %s and must be rebuilt.\n"
                "Build one with: python tools/build_reference.py"
                % (bundle, ', '.join(stale)))
        meta = dict(self._db.execute('SELECT key, value FROM meta').fetchall())
        self.release = meta.get('release', 'unknown')
        self.assembly = meta.get('assembly', 'GRCh38')
        self._genome = None

    def describe(self):
        return 'local bundle %s (Ensembl release %s)' % (self.bundle, self.release)

    def _fasta(self):
        if self._genome is None:
            from pyfaidx import Fasta
            # sequence_always_upper guards against a bundle built from a
            # soft-masked assembly, where lowercase bases would pass the
            # script's sequence check but match no 'C'/'A' edit.
            self._genome = Fasta(self.genome_path, as_raw=True,
                                 sequence_always_upper=True)
        return self._genome

    def _record(self, tr):
        row = self._db.execute(
            'SELECT * FROM transcript WHERE transcript_id = ?', (tr,)).fetchone()
        if row is None:
            raise TranscriptNotFound(tr, 'not in bundle %s' % os.path.basename(self.bundle))
        return row

    def lookup(self, tr):
        row = self._record(tr)
        return {
            'id': tr,
            'display_name': row['display_name'],
            'assembly_name': self.assembly,
            'strand': row['strand'],
            'seq_region_name': row['seq_region'],
            'Parent': row['gene_id'],
            'biotype': row['biotype'],
            'start': row['start'],
            'end': row['end'],
            'Exon': [{'start': s, 'end': e} for s, e in json.loads(row['exons'])],
            # Beyond what REST returns: which transcript a gene search should
            # preselect. __init__ has already refused a bundle without these.
            'gene_name': row['gene_name'],
            'mane_select': bool(row['mane_select']),
            'ensembl_canonical': bool(row['ensembl_canonical']),
        }

    def cds_mappings(self, tr, length):
        row = self._record(tr)
        if not row['cds']:
            return None
        # /map/cds/{tr}/1..length stops after `length` coding bases. Callers pass
        # the genomic span, which for a single-exon transcript with no UTR can
        # fall one base short of the CDS, so honour the limit.
        mappings = []
        used = 0
        for start, end in json.loads(row['cds']):
            span = end - start + 1
            if used + span > length:
                keep = length - used
                if keep <= 0:
                    break
                # trim from the 3' end, which is the high coordinate on + strand
                if row['strand'] == 1:
                    end = start + keep - 1
                else:
                    start = end - keep + 1
                span = keep
            mappings.append({'start': start, 'end': end, 'strand': row['strand']})
            used += span
            if used >= length:
                break
        return mappings or None

    def genomic_sequence(self, tr, flank=40):
        row = self._record(tr)
        genome = self._fasta()
        contig = row['seq_region']
        if contig not in genome:
            return ''
        # Ensembl clamps the flank at a contig edge rather than failing.
        start = max(1, row['start'] - flank)
        end = min(len(genome[contig]), row['end'] + flank)
        seq = genome[contig][start - 1:end]
        if row['strand'] == -1:
            seq = reverse_complement(seq)
        return seq

    def protein_sequence(self, tr):
        return self._record(tr)['protein_sequence'] or ''

    def cds_sequence(self, tr):
        return self._record(tr)['cds_sequence'] or ''


def empty_variants():
    """The no-variants case, shaped like `variants_for_gene` output.

    Used for nucleotide input and for --no-clinvar. `ClinVar_SNP_Position` has to be
    an integer column so `get_snps` can match it against genomic positions.
    """
    frame = pd.DataFrame(columns=list(VARIANT_COLUMNS), dtype=object)
    return frame.astype({'ClinVar_SNP_Position': 'int64'})


class ClinVarSource(object):
    """ClinVar variants on disk, from a database built by tools/build_clinvar.py.

    A design run only ever asks for one gene's variants at a time, so this queries
    an indexed table rather than loading the 4.1-million-row export.
    """

    def __init__(self, db_path):
        if not os.path.exists(db_path):
            raise ClinVarNotFound(
                "ClinVar database '%s' does not exist.\n"
                "Build one with: python tools/build_clinvar.py" % db_path)
        self.db_path = db_path
        self._db = sqlite3.connect('file:%s?mode=ro' % db_path, uri=True)
        meta = dict(self._db.execute('SELECT key, value FROM meta').fetchall())
        self.rows = meta.get('rows', 'unknown')
        self._cache = {}

    def describe(self):
        """One line naming the variant data, recorded in the run README."""
        return 'ClinVar database %s (%s variants)' % (self.db_path, self.rows)

    def variants_for_gene(self, symbol):
        """The gene's variants, in the shape design_sgrnas expects."""
        if symbol not in self._cache:
            if len(self._cache) >= VARIANT_CACHE_SIZE:
                # dicts keep insertion order, so this drops the oldest symbol
                self._cache.pop(next(iter(self._cache)))
            self._cache[symbol] = pd.read_sql(
                'SELECT * FROM %s WHERE GeneSymbol = ?' % VARIANT_TABLE,
                self._db, params=(symbol,))
        return self._cache[symbol]


def find_clinvar_db(clinvar_db, refdata='refdata'):
    """Resolves --clinvar-db to a database file.

    Accepts the file itself, or falls back to the newest clinvar-<date>.db under
    --refdata. The names carry ISO dates, so they sort by age.
    """
    if clinvar_db:
        return clinvar_db
    candidates = sorted(glob.glob(os.path.join(refdata, CLINVAR_DB_GLOB)))
    if not candidates:
        raise ClinVarNotFound(
            "No ClinVar database under '%s'.\n"
            "Build one with: python tools/build_clinvar.py\n"
            "Or skip ClinVar annotation with --no-clinvar." % refdata)
    return candidates[-1]


def find_bundle(refdata):
    """Resolves --refdata to a bundle directory.

    Accepts either the bundle itself or a parent holding ensembl-<release>
    directories, in which case the highest release wins.
    """
    if not missing_bundle_files(refdata):
        return refdata
    candidates = sorted(
        (d for d in glob.glob(os.path.join(refdata, 'ensembl-*'))
         if not missing_bundle_files(d)),
        key=lambda d: _release_key(os.path.basename(d)))
    if not candidates:
        raise BundleNotFound(
            "No reference bundle under '%s'.\n"
            "Build one with: python tools/build_reference.py" % refdata)
    return candidates[-1]


def _release_key(name):
    m = re.search(r'(\d+)$', name)
    return int(m.group(1)) if m else -1


def get_source(name='rest', refdata='refdata'):
    """Builds the source named on the command line."""
    if name == 'rest':
        return EnsemblRestSource()
    if name == 'local':
        return LocalSource(find_bundle(refdata))
    raise ValueError("Unknown transcript source '%s'" % name)
