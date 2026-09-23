"""Running design jobs in a process pool.

Design is CPU-bound, so it runs in `ProcessPoolExecutor` workers. Each worker opens
its own `LocalSource` and `ClinVarSource` in the initializer: both hold sqlite3
handles, which cannot be pickled across the process boundary, and both are cheap.

The worker writes the result objects itself rather than returning rows. TTN is 25,662
guides, and sending that back through the pool's pickle channel would make the parent
pay for the size of the answer.
"""
import csv
import dataclasses
import gzip
import io
import json
import logging
import multiprocessing
import signal
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool

from bedesign import (ANNOTATION_COLUMNS, DESIGN_COLUMNS, ENGINE_VERSION,
                      ERROR_COLUMNS, design_transcript)
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
        'params': dataclasses.asdict(params),
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


class ClientBusy(Exception):
    """This client already has a design running; it should wait for that one."""


class JobPool(object):
    """Tracks which keys are running, which clients started them, and which failed.

    A single lock covers every table. Checking and inserting under it is what makes
    two simultaneous requests for one uncached key run one job rather than two.

    One running job per client. The rate limiter counts job starts, but a design
    can hold a worker for the whole job timeout, so without this one client could
    hold every worker while staying inside its rate limit.

    The executor is replaced when a worker dies. ProcessPoolExecutor does not
    recover from a worker killed outright (the OOM killer, a crash in native code):
    every later submit raises BrokenProcessPool, forever.
    """

    def __init__(self, references, settings, after_job=None):
        self.settings = settings
        # Called with no arguments after each job that wrote a result, outside the
        # lock. The app uses it to keep the results directory under its budget.
        self._after_job = after_job
        self._initargs = (references.bundle, references.clinvar_db, settings.results_dir)
        self._lock = threading.Lock()
        self._running = {}
        self._clients = {}
        self._failed = {}
        # Set when a worker dies, cleared when a job next finishes without one
        # dying. What /healthz reports.
        self._crashed = False
        self._pool = self._new_pool()

    def _new_pool(self):
        return ProcessPoolExecutor(
            max_workers=self.settings.max_jobs,
            # 'spawn' explicitly, rather than whatever the platform defaults to.
            # The parent is a threaded ASGI server, and forking a process that
            # holds locks is how a worker deadlocks before it runs anything.
            mp_context=multiprocessing.get_context('spawn'),
            initializer=init_worker,
            initargs=self._initargs)

    def _replace_pool(self):
        """Called with the lock held, only for a pool that is already broken.

        The old pool is dropped, not shut down. A broken executor has already
        terminated its workers and stopped its manager thread. It also fails its
        pending futures while holding its own non-reentrant shutdown lock, so
        calling shutdown() from their done-callbacks -- which is where `_finished`
        runs -- deadlocks, and with it every request waiting on this pool's lock.
        """
        self._pool = self._new_pool()

    def healthy(self):
        """False from a worker dying until a job next finishes normally."""
        with self._lock:
            return not self._crashed

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

    def client_busy(self, client):
        """Whether this client has a job running. Lets the caller turn a second
        job away before charging the client a rate-limit token for it."""
        with self._lock:
            return client in self._clients

    def submit(self, key, transcript_id, params, client=None):
        """Starts a job unless one is already running for this key.

        Returns True if this call started it, False if it joined one already in
        flight -- joining is free, whoever started the job. Raises ClientBusy when
        `client` already has a different job running, and PoolBusy when every
        worker is occupied.
        """
        with self._lock:
            if key in self._running:
                return False
            if client is not None and client in self._clients:
                raise ClientBusy()
            if len(self._running) >= self.settings.max_jobs:
                # No queue: the design doc's model is that an unserved request
                # retries, so a backlog never builds up behind a slow gene.
                raise PoolBusy()
            self._failed.pop(key, None)
            args = (run_job, key, transcript_id, params, self.settings.job_timeout)
            try:
                future = self._pool.submit(*args)
            except BrokenProcessPool:
                # A worker died while no job was running, so no future reported
                # it. The replacement is what this request runs on.
                log.error('design pool was broken; replacing it')
                self._replace_pool()
                future = self._pool.submit(*args)
            pool = self._pool
            self._running[key] = future
            if client is not None:
                self._clients[client] = key
        future.add_done_callback(lambda f: self._finished(key, client, pool, f))
        return True

    def _finished(self, key, client, pool, future):
        message = None
        crashed = False
        try:
            future.result()
        except JobTimeout:
            message = ('This design took longer than %d seconds and was stopped.'
                       % self.settings.job_timeout)
        except BrokenProcessPool:
            log.error('a design worker died running key %s; replacing the pool', key)
            message = 'This design failed. Reload to try again.'
            crashed = True
        except Exception:
            # The detail goes to the log; the page gets a fixed string.
            log.exception('design job failed for key %s', key)
            message = 'This design failed. Reload to try again.'
        with self._lock:
            self._running.pop(key, None)
            if client is not None and self._clients.get(client) == key:
                del self._clients[client]
            if message is not None:
                self._failed[key] = (message, time.time())
            self._crashed = crashed
            # Every job on a broken pool fails at once. Only the first to get
            # here replaces it; the rest find a newer pool already in place.
            if crashed and pool is self._pool:
                self._replace_pool()
        if message is None and self._after_job is not None:
            try:
                self._after_job()
            except Exception:
                # Housekeeping; a failure here must not look like a failed design.
                log.exception('after-job hook failed for key %s', key)

    def shutdown(self):
        self._pool.shutdown(wait=False, cancel_futures=True)
