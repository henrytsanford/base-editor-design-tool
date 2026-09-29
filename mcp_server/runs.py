"""Designs, run once and addressed by what produced them.

A run is stored the way the web service stores one -- three gzipped TSVs and a
manifest written last, under service.storage -- so the same bytes back a download, a
table and a tool result. The handle is service.cachekey's digest over the resolved
parameters, which makes an identical second request free and means a handle can never
name a result computed under parameters the caller has forgotten.

The web service writes into this same store, under this same key, so a transcript
designed there and designed here are one object. Its manifest is the smaller, older
shape; `_manifest` fills in what it does not carry rather than failing on a key while
the result sits on disk.
"""
import contextlib
import csv
import dataclasses
import gzip
import hashlib
import io
import json
import logging
import threading
import time

from bedesign import (ANNOTATION_COLUMNS, ANNOTATIONS_FILE, DESIGN_COLUMNS,
                      DESIGNS_FILE, ENGINE_VERSION, ERROR_COLUMNS, ERRORS_FILE,
                      design_sequence, design_transcript, tsv)
from bedesign.transcript_source import ClinVarSource, LocalSource

from service.cachekey import (MANIFEST, RESULT_FILES, RESULTS, cache_key,
                              manifest_key, result_key)
from service.results import ResultCache
from service.storage import LocalStorage
from service.tables import load_table

log = logging.getLogger(__name__)

# The three outputs a run writes: the columns each holds, and the name the CLI would
# have given the file. One table, because the same three names are otherwise spelled
# out wherever a file is written, validated, exported or named.
OUTPUTS = {
    'designs': (DESIGN_COLUMNS, DESIGNS_FILE),
    'errors': (ERROR_COLUMNS, ERRORS_FILE),
    'clinvar': (ANNOTATION_COLUMNS, ANNOTATIONS_FILE),
}
# The storage layer owns the same vocabulary; a name here that it does not know would
# be a file nothing could ever read back.
assert set(OUTPUTS) == set(RESULT_FILES), 'output names disagree with storage'


def _file_stats(raw):
    """A written file's size, line count and digest.

    Measured over the text form -- the bytes `export_run` hands back -- rather than
    the stored gzip, so the figures describe the same content whichever encoding a
    caller asks for.
    """
    text = tsv.as_text(raw)
    encoded = text.encode('utf-8')
    return {'bytes': len(encoded), 'lines': text.count('\n'),
            'sha256': hashlib.sha256(encoded).hexdigest()}


def filename_for(name, label):
    """What the CLI would have called this output for this run."""
    template = OUTPUTS[name][1]
    return template % label if '%s' in template else template


class RunNotFound(KeyError):
    """A handle naming a result that is not stored any more."""


class ResultTooLarge(ValueError):
    """A result with more guides than the table layer will parse."""


@dataclasses.dataclass(frozen=True)
class Target:
    """What a run designs over: a transcript, or a sequence given inline."""
    kind: str
    transcript_id: str = ''
    name: str = ''
    sequence: str = ''

    @property
    def label(self):
        """What the output files are named after, as the CLI names them."""
        return self.transcript_id if self.kind == 'transcript' else self.name

    @property
    def key(self):
        """The part of the cache key that says what was designed.

        A sequence is keyed by its digest rather than by itself: the key is canonical
        JSON that gets hashed anyway, and a 100 kb sequence has no business being
        built into a string first.
        """
        if self.kind == 'transcript':
            return self.transcript_id
        digest = hashlib.sha256(
            ('%s\n%s' % (self.name, self.sequence)).encode('utf-8')).hexdigest()
        return 'seq:%s' % digest

    @property
    def files(self):
        """Which of the three files this kind of run writes.

        Nucleotide input has no gene to look up, so it never annotates -- the CLI does
        not create the file at all, and neither does this.
        """
        return ('designs', 'errors', 'clinvar') if self.kind == 'transcript' \
            else ('designs', 'errors')


def _design(target, source, clinvar, params):
    """The only call into the engine: the transcript / sequence fork, in one place.

    Nothing here redirects stdout. The engine and the transcript source log their
    diagnostics rather than printing them, which is what keeps them off the stdio
    transport's stream -- a redirect could not, since sys.stdout is process-global
    and this server answers several calls at once.
    """
    if target.kind == 'transcript':
        return design_transcript(source, clinvar, target.transcript_id, params)
    return design_sequence(target.name, target.sequence, params)


class RunStore(object):
    """Runs designs, stores them, and hands back parsed tables."""

    def __init__(self, settings, references):
        self.settings = settings
        self.references = references
        self.storage = LocalStorage(settings.results_dir)
        self.tables = ResultCache(settings.table_cache_rows,
                                  settings.table_cache_frames)
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
        existing = self._manifest(run_id)
        if existing is not None:
            self.storage.touch(manifest_key(run_id))
            return run_id, existing

        started = time.time()
        designs, errors, annotations = _design(
            target, self.source, self.clinvar if target.kind == 'transcript' else None,
            params)
        rows = {'designs': designs, 'errors': errors, 'clinvar': annotations}
        stats = {}
        for name in target.files:
            raw = tsv.tsv_bytes(OUTPUTS[name][0], rows[name])
            self.storage.put(result_key(run_id, name), tsv.gzipped(raw))
            # Measured here because the uncompressed bytes are in hand exactly once.
            # Recording them lets export_run answer for a file's size, line count and
            # digest without decompressing it -- including when it has to refuse an
            # oversize export, which otherwise pays for the whole file to say no.
            stats[name] = _file_stats(raw)

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

    def _manifest(self, run_id):
        try:
            stored = json.loads(self.storage.get(manifest_key(run_id)))
        except (KeyError, ValueError):
            return None
        return self._filled(run_id, stored)

    def _filled(self, run_id, manifest):
        """A stored manifest with the fields only this front door writes.

        The web service designs into the same store under the same key, so a run it
        computed is a run this can serve -- but its manifest predates these fields.
        Filling them in on read means a result that is sitting on disk is returned,
        rather than raising a KeyError about a shape difference the caller cannot see
        and did not cause.
        """
        if 'run_id' in manifest and 'files' in manifest:
            return manifest
        transcript_id = manifest.get('transcript_id', '')
        filled = dict(manifest)
        filled.setdefault('run_id', run_id)
        filled.setdefault('kind', 'transcript' if transcript_id else 'sequence')
        filled.setdefault('label', transcript_id)
        filled.setdefault('gene', self.references.gene_name(transcript_id)
                          if transcript_id else '')
        filled.setdefault('stats', {})
        filled.setdefault('files', [name for name in OUTPUTS
                                    if self.storage.exists(result_key(run_id, name))])
        return filled

    def manifest(self, run_id):
        """The manifest for a handle, or RunNotFound.

        A handle whose result was swept is not an error in the caller's reasoning, so
        the message says what to do about it rather than what went wrong.
        """
        manifest = self._manifest(run_id)
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
        return recorded or _file_stats(gzip.decompress(self.stored(manifest, name)))

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
