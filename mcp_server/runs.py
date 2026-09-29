"""Designs, run once and addressed by what produced them.

A run is stored the way the web service stores one -- three gzipped TSVs and a
manifest written last, under service.storage -- so the same bytes back a download, a
table and a tool result. The handle is service.cachekey's digest over the resolved
parameters, which makes an identical second request free and means a handle can never
name a result computed under parameters the caller has forgotten.

The web service writes into this same store, under this same key, so a transcript
designed there and designed here are one object. What each writer records differs a
little, so `_filled` supplies the rest on read rather than failing on a key while the
result sits on disk.
"""
import contextlib
import csv
import dataclasses
import gzip
import io
import json
import logging
import threading
import time

from bedesign import (ANNOTATIONS_FILE, DESIGNS_FILE, ENGINE_VERSION, ERRORS_FILE,
                      tsv)
from bedesign.transcript_source import ClinVarSource, LocalSource

from service.cachekey import (MANIFEST, RESULT_FILES, RESULTS, cache_key,
                              manifest_key, result_key)
from service.jobs import COLUMNS, design_target, file_stats
from service.results import ResultCache
from service.storage import LocalStorage
from service.tables import load_table

log = logging.getLogger(__name__)

# What the CLI would have called each output. The columns each holds live beside the
# worker that writes them, in service.jobs.COLUMNS.
FILENAMES = {'designs': DESIGNS_FILE, 'errors': ERRORS_FILE,
             'clinvar': ANNOTATIONS_FILE}
# The storage layer owns the same vocabulary; a name here that it does not know would
# be a file nothing could ever read back.
assert set(FILENAMES) == set(RESULT_FILES), 'output names disagree with storage'


def filename_for(name, label):
    """What the CLI would have called this output for this run."""
    template = FILENAMES[name]
    return template % label if '%s' in template else template


class RunNotFound(KeyError):
    """A handle naming a result that is not stored any more."""


class ResultTooLarge(ValueError):
    """A result with more guides than the table layer will parse."""


class ServiceBusy(Exception):
    """The service could not start or finish this design now; retrying is the answer."""


class DesignFailed(Exception):
    """A design ran and produced nothing."""


class RunStore(object):
    """Runs designs, stores them, and hands back parsed tables."""

    def __init__(self, settings, references, storage=None, tables=None):
        """`storage` and `tables` are the web app's own when this is mounted on it.

        Sharing them is not an optimisation. Two LocalStorage objects over one
        directory hold two independent evict locks, so the two would sweep the results
        tree against each other; two ResultCaches would each hold the full
        `table_cache_rows` budget, doubling the memory that setting exists to bound.
        Shared, a result parsed for a browser table is the same frame a query_guides
        page reads. The stdio server passes neither and gets its own, as before.

        `is not None` rather than `or`: an empty ResultCache is truthy today, but what
        is being asked is whether one was given, and that should not depend on a
        __bool__ nobody wrote.
        """
        self.settings = settings
        self.references = references
        self.storage = (storage if storage is not None
                        else LocalStorage(settings.results_dir))
        self.tables = (tables if tables is not None
                       else ResultCache(settings.table_cache_rows,
                                        settings.table_cache_frames))
        self._local = threading.local()

    # Sources are per-thread, for the reason service/references.py gives at length:
    # a sqlite3 connection belongs to the thread that opened it, and the MCP server
    # answers each call on a worker thread from a pool. One shared LocalSource raises
    # ProgrammingError on the second call to arrive on a new thread. Opening a
    # read-only connection costs about 0.2 ms and the pool is bounded, so one per
    # thread is cheaper than serialising every design behind a lock.

    @property
    def source(self):
        """This thread's transcript source, opened on first use."""
        source = getattr(self._local, 'source', None)
        if source is None:
            source = self._local.source = LocalSource(self.references.bundle)
        return source

    @property
    def clinvar(self):
        """This thread's ClinVar database, opened on first use."""
        clinvar = getattr(self._local, 'clinvar', None)
        if clinvar is None:
            clinvar = self._local.clinvar = ClinVarSource(self.references.clinvar_db)
        return clinvar

    # ---- running -------------------------------------------------------------

    def design(self, target, params):
        """Designs `target`, or returns the manifest of an identical earlier run.

        Content-addressed, so asking the same question twice in one conversation
        costs one design. The manifest is written last, which is what makes a handle
        that resolves a handle whose files are all there.
        """
        run_id = cache_key(target.key, params, self.references.release,
                           self.references.clinvar_version, ENGINE_VERSION)
        existing = self.lookup(run_id)
        if existing is not None:
            return run_id, existing

        started = time.time()
        designs, errors, annotations = design_target(
            target, self.source, self.clinvar if target.kind == 'transcript' else None,
            params)
        rows = {'designs': designs, 'errors': errors, 'clinvar': annotations}
        stats = {}
        for name in target.files:
            raw = tsv.tsv_bytes(COLUMNS[name], rows[name])
            self.storage.put(result_key(run_id, name), tsv.gzipped(raw))
            # Measured here because the uncompressed bytes are in hand exactly once.
            # Recording them lets export_run answer for a file's size, line count and
            # digest without decompressing it -- including when it has to refuse an
            # oversize export, which otherwise pays for the whole file to say no.
            stats[name] = file_stats(raw)

        manifest = {
            'run_id': run_id,
            'kind': target.kind,
            'transcript_id': target.transcript_id,
            'gene': (self.references.gene_name(target.transcript_id)
                     if target.kind == 'transcript' else ''),
            'label': target.label,
            'params': dataclasses.asdict(params),
            'engine_version': ENGINE_VERSION,
            'reference': self.source.describe(),
            'clinvar': (self.clinvar.describe()
                        if target.kind == 'transcript' else 'none'),
            'designs': len(designs),
            'errors': len(errors),
            'annotations': len(annotations) if target.kind == 'transcript' else 0,
            'files': list(target.files),
            'stats': stats,
            'runtime_seconds': round(time.time() - started, 3),
        }
        # Written last: a manifest present means every file beside it is complete.
        self.storage.put(manifest_key(run_id),
                         json.dumps(manifest).encode('utf-8'))
        self.storage.evict(RESULTS, MANIFEST, self.settings.results_max_mb * 2 ** 20)
        return run_id, manifest

    # ---- reading back --------------------------------------------------------

    def lookup(self, run_id):
        """The manifest for a handle, or None if nothing complete is stored.

        The manifest is written last, so its presence is what says the result behind
        it is whole. A hit is touched, because a view is a use: it is what keeps a
        result a caller is still reading from being the next one evicted.
        """
        try:
            stored = json.loads(self.storage.get(manifest_key(run_id)))
        except (KeyError, ValueError):
            return None
        self.storage.touch(manifest_key(run_id))
        return self._filled(run_id, stored)

    def _filled(self, run_id, manifest):
        """A stored manifest with whatever its writer did not record.

        Two writers put manifests in this store. The pool worker records everything it
        can measure, but resolving a gene symbol there would cost a SELECT * over a row
        holding the CDS and protein sequences, so the symbol is filled in here off the
        index. An older manifest that predates the other fields is filled the same way,
        which is what lets a result sitting on disk be returned rather than raising
        about a shape difference the caller cannot see and did not cause.
        """
        transcript_id = manifest.get('transcript_id', '')
        filled = dict(manifest)
        filled.setdefault('run_id', run_id)
        filled.setdefault('kind', 'transcript' if transcript_id else 'sequence')
        filled.setdefault('label', transcript_id)
        if 'gene' not in filled:
            filled['gene'] = (self.references.gene_name(transcript_id)
                              if transcript_id else '')
        filled.setdefault('stats', {})
        filled.setdefault('files', [name for name in COLUMNS
                                    if self.storage.exists(result_key(run_id, name))])
        return filled

    def manifest(self, run_id):
        """The manifest for a handle, or RunNotFound.

        A handle whose result was swept is not an error in the caller's reasoning, so
        the message says what to do about it rather than what went wrong.
        """
        manifest = self.lookup(run_id)
        if manifest is None:
            raise RunNotFound(
                'No stored result for that run_id. It may have been cleared; '
                'run design_guides again to recompute it.')
        return manifest

    def stored(self, manifest, name):
        """One stored file's compressed bytes, checked against what the run wrote."""
        if name not in manifest.get('files', ()):
            raise RunNotFound(
                'This run has no %s file. A sequence run has no gene to look up, '
                'so it produces no ClinVar annotations.' % name)
        return self.storage.get(result_key(manifest['run_id'], name))

    def stats(self, manifest, name):
        """A file's size, line count and digest, measured when it was written."""
        recorded = manifest.get('stats', {}).get(name)
        # A run the web service wrote carries no measurements, so they are taken now.
        return recorded or file_stats(gzip.decompress(self.stored(manifest, name)))

    def table(self, run_id):
        """The parsed designs table for a handle, with coverage attached.

        Bounded the way the web app bounds it: past `table_max_rows` a frame is large
        enough to matter, and the answer is the whole file rather than a page of it.
        """
        manifest = self.manifest(run_id)
        if manifest['designs'] > self.settings.table_max_rows:
            raise ResultTooLarge(
                '%d guides is too many to page through; use export_run to take the '
                'whole file, or design with a more specific PAM.'
                % manifest['designs'])

        transcript_id = (manifest['transcript_id']
                         if manifest['kind'] == 'transcript' else '')
        self.storage.touch(manifest_key(run_id))
        return self.tables.get(run_id, lambda: load_table(
            self.storage, self.references, run_id, transcript_id,
            manifest['params']['intron_buffer']))

    @contextlib.contextmanager
    def annotation_reader(self, manifest):
        """The annotation rows as a stream, header already consumed.

        Streamed rather than materialised: the file runs to several megabytes on a
        large gene and a page is fifty rows, so building the whole list to slice it
        would cost more memory than the design did.

        Read with csv rather than pandas, for the reason CLAUDE.md gives: a row that
        matched no ClinVar SNP is truncated after 'Mutation category', and read_csv
        would pad every one of them out to 28 fields with blanks the file never had.
        """
        blob = self.stored(manifest, 'clinvar')
        # newline=None, so the CRLF the file carries is translated the same way
        # tsv.as_text translates it and a row's last field keeps no stray return.
        with gzip.open(io.BytesIO(blob), 'rt', encoding='utf-8') as fh:
            reader = csv.reader(fh, delimiter='\t')
            next(reader, None)
            yield reader
