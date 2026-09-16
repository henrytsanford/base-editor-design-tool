"""Running design jobs in a process pool.

Design is CPU-bound, so it runs in `ProcessPoolExecutor` workers. Each worker opens
its own `LocalSource` and `ClinVarSource` in the initializer: both hold sqlite3
handles, which cannot be pickled across the process boundary, and both are cheap.

The worker writes the result objects itself rather than returning rows. TTN is 25,662
guides, and sending that back through the pool's pickle channel would make the parent
pay for the size of the answer.
"""
import csv
import gzip
import io
import json
import logging
import multiprocessing
import signal
import threading
import time
from concurrent.futures import ProcessPoolExecutor

from bedesign import (ANNOTATION_COLUMNS, DESIGN_COLUMNS, ENGINE_VERSION,
                      ERROR_COLUMNS, DesignParams, design_transcript)
from bedesign.transcript_source import ClinVarSource, LocalSource

from .cachekey import manifest_key, result_key
from .storage import LocalStorage

log = logging.getLogger(__name__)

# Worker-process state, built once per process by `init_worker`.
_source = None
_clinvar = None
_storage = None


class JobTimeout(Exception):
    """A design job passed its wall-clock cap."""


def init_worker(bundle, clinvar_db, results_dir):
    global _source, _clinvar, _storage
    _source = LocalSource(bundle)
    _clinvar = ClinVarSource(clinvar_db)
    _storage = LocalStorage(results_dir)


def _tsv_gz(columns, rows):
    """A gzipped TSV, byte-identical for identical rows.

    mtime=0 for the same reason the golden fixtures use it: the bytes should depend
    on the designs, not on when they were computed.
    """
    raw = io.BytesIO()
    with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as gz:
        text = io.TextIOWrapper(gz, encoding='utf-8', newline='')
        writer = csv.writer(text, delimiter='\t')
        writer.writerow(columns)
        writer.writerows(rows)
        text.flush()
        text.detach()
    return raw.getvalue()


def run_job(key, transcript_id, params, timeout):
    """Designs one transcript and writes its results. Runs in a worker process.

    The wall-clock cap is enforced here rather than in the parent because
    ProcessPoolExecutor cannot cancel a future that has already started: a
    parent-side timeout would return to the user while the worker kept a core busy.
    """
    started = time.time()

    def expired(signum, frame):
        raise JobTimeout('design exceeded %d seconds' % timeout)

    previous = signal.signal(signal.SIGALRM, expired)
    signal.alarm(timeout)
    try:
        designs, errors, annotations = design_transcript(
            _source, _clinvar, transcript_id, params)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)

    runtime = time.time() - started
    _storage.put(result_key(key, 'designs'), _tsv_gz(DESIGN_COLUMNS, designs))
    _storage.put(result_key(key, 'errors'), _tsv_gz(ERROR_COLUMNS, errors))
    _storage.put(result_key(key, 'clinvar'), _tsv_gz(ANNOTATION_COLUMNS, annotations))
    manifest = {
        'transcript_id': transcript_id,
        'params': {f: getattr(params, f) for f in
                   ('pam', 'window', 'sg_len', 'edit', 'intron_buffer', 'filter_gc')},
        'engine_version': ENGINE_VERSION,
        'reference': _source.describe(),
        'clinvar': _clinvar.describe(),
        'designs': len(designs),
        'errors': len(errors),
        'annotations': len(annotations),
        'runtime_seconds': round(runtime, 3),
    }
    # Written last, so its presence means the result behind it is complete.
    _storage.put(manifest_key(key),
                 json.dumps(manifest, sort_keys=True).encode('utf-8'))
    return {'designs': len(designs), 'runtime_seconds': round(runtime, 3)}


class PoolBusy(Exception):
    """Every worker is occupied; the caller should retry rather than queue."""


class JobPool(object):
    """Tracks which keys are running, and which failed recently.

    A single lock covers both dicts. Checking and inserting under it is what makes
    two simultaneous requests for one uncached key run one job rather than two.
    """

    def __init__(self, references, settings):
        self.settings = settings
        self._lock = threading.Lock()
        self._running = {}
        self._failed = {}
        self._pool = ProcessPoolExecutor(
            max_workers=settings.max_jobs,
            # 'spawn' explicitly, rather than whatever the platform defaults to.
            # The parent is a threaded ASGI server, and forking a process that
            # holds locks is how a worker deadlocks before it runs anything.
            mp_context=multiprocessing.get_context('spawn'),
            initializer=init_worker,
            initargs=(references.bundle, references.clinvar_db, settings.results_dir))

    def failure(self, key):
        """The recent failure for a key, if it has not aged out yet."""
        with self._lock:
            entry = self._failed.get(key)
            if entry is None:
                return None
            message, when = entry
            if time.time() - when > self.settings.failure_ttl:
                del self._failed[key]
                return None
            return message

    def running(self, key):
        with self._lock:
            return key in self._running

    def submit(self, key, transcript_id, params):
        """Starts a job unless one is already running for this key.

        Returns True if this call started it, False if it joined one already in
        flight. Raises PoolBusy when every worker is occupied.
        """
        with self._lock:
            if key in self._running:
                return False
            if len(self._running) >= self.settings.max_jobs:
                # No queue: the design doc's model is that an unserved request
                # retries, so a backlog never builds up behind a slow gene.
                raise PoolBusy()
            self._failed.pop(key, None)
            future = self._pool.submit(
                run_job, key, transcript_id, params, self.settings.job_timeout)
            self._running[key] = future
        future.add_done_callback(lambda f: self._finished(key, f))
        return True

    def _finished(self, key, future):
        message = None
        try:
            future.result()
        except JobTimeout:
            message = ('This design took longer than %d seconds and was stopped.'
                       % self.settings.job_timeout)
        except Exception:
            # The detail goes to the log; the page gets a fixed string.
            log.exception('design job failed for key %s', key)
            message = 'This design failed. Reload to try again.'
        with self._lock:
            self._running.pop(key, None)
            if message is not None:
                self._failed[key] = (message, time.time())

    def shutdown(self):
        self._pool.shutdown(wait=False, cancel_futures=True)
