"""The job table: deduplication, the busy path, failures and the timeout.

The pool is replaced with a fake so these run offline and deterministically. What is
under test is the bookkeeping in JobPool, not the design engine.
"""
import time
from concurrent.futures import Future

import pytest

import service.jobs as jobs
from service.config import Settings
from service.jobs import JobPool, JobTimeout, PoolBusy


class FakeExecutor(object):
    def __init__(self, **kwargs):
        self.calls = []

    def submit(self, fn, *args):
        future = Future()
        self.calls.append((args, future))
        return future

    def shutdown(self, **kwargs):
        pass


class FakeReferences(object):
    bundle = 'bundle'
    clinvar_db = 'clinvar.db'


@pytest.fixture
def pool(monkeypatch):
    monkeypatch.setattr(jobs, 'ProcessPoolExecutor',
                        lambda **kwargs: FakeExecutor(**kwargs))
    return JobPool(FakeReferences(), Settings(max_jobs=2, failure_ttl=300))


def test_two_requests_for_one_key_start_one_job(pool):
    """The milestone's own criterion: requesting an uncached result twice runs it once."""
    assert pool.submit('key', 'ENST1', None) is True
    assert pool.submit('key', 'ENST1', None) is False
    assert len(pool._pool.calls) == 1


def test_a_running_key_is_reported_as_running(pool):
    pool.submit('key', 'ENST1', None)
    assert pool.running('key') is True
    assert pool.running('other') is False


def test_a_full_pool_refuses_rather_than_queueing(pool):
    pool.submit('a', 'ENST1', None)
    pool.submit('b', 'ENST2', None)
    with pytest.raises(PoolBusy):
        pool.submit('c', 'ENST3', None)


def test_a_finished_job_leaves_the_running_table(pool):
    pool.submit('key', 'ENST1', None)
    _, future = pool._pool.calls[0]
    future.set_result({'designs': 1})
    assert pool.running('key') is False
    assert pool.failure('key') is None


def test_a_finished_job_frees_its_slot(pool):
    pool.submit('a', 'ENST1', None)
    pool.submit('b', 'ENST2', None)
    pool._pool.calls[0][1].set_result({})
    pool.submit('c', 'ENST3', None)  # would raise PoolBusy if the slot were held


def test_a_failure_is_remembered_without_leaking_detail(pool):
    pool.submit('key', 'ENST1', None)
    pool._pool.calls[0][1].set_exception(RuntimeError('/refdata/secret.db is missing'))
    message = pool.failure('key')
    assert message == 'This design failed. Reload to try again.'
    assert 'secret' not in message and 'refdata' not in message


def test_a_timeout_says_so(pool):
    pool.submit('key', 'ENST1', None)
    pool._pool.calls[0][1].set_exception(JobTimeout('too slow'))
    assert 'longer than' in pool.failure('key')


def test_a_failure_is_forgotten_so_a_reload_retries(monkeypatch):
    monkeypatch.setattr(jobs, 'ProcessPoolExecutor',
                        lambda **kwargs: FakeExecutor(**kwargs))
    pool = JobPool(FakeReferences(), Settings(max_jobs=2, failure_ttl=0))
    pool.submit('key', 'ENST1', None)
    pool._pool.calls[0][1].set_exception(RuntimeError('boom'))
    time.sleep(0.01)
    assert pool.failure('key') is None


def test_resubmitting_clears_a_previous_failure(pool):
    pool.submit('key', 'ENST1', None)
    pool._pool.calls[0][1].set_exception(RuntimeError('boom'))
    assert pool.failure('key') is not None
    pool.submit('key', 'ENST1', None)
    assert pool.failure('key') is None


def test_the_worker_timeout_interrupts_a_slow_design(monkeypatch):
    """The cap is enforced inside the worker: ProcessPoolExecutor cannot cancel a
    future that has already started, so a parent-side timeout would leave a core busy."""
    def slow(*args, **kwargs):
        time.sleep(30)

    monkeypatch.setattr(jobs, 'design_transcript', slow)
    started = time.time()
    with pytest.raises(JobTimeout):
        jobs.run_job('key', 'ENST1', None, timeout=1)
    assert time.time() - started < 5


def test_the_alarm_is_cleared_when_a_design_finishes(monkeypatch):
    """A leftover alarm would fire during an unrelated later request."""
    import signal
    monkeypatch.setattr(jobs, 'design_transcript', lambda *a, **k: ([], [], []))
    monkeypatch.setattr(jobs, '_storage', type('S', (), {'put': lambda *a: None})())
    monkeypatch.setattr(jobs, '_source', type('S', (), {'describe': lambda s: 'x'})())
    monkeypatch.setattr(jobs, '_clinvar', type('C', (), {'describe': lambda s: 'y'})())
    jobs.run_job('key', 'ENST1', __import__('bedesign').DesignParams(), timeout=60)
    assert signal.alarm(0) == 0


def test_tsv_bytes_do_not_depend_on_when_they_were_written():
    """mtime=0, for the reason the golden fixtures use it."""
    first = jobs._tsv_gz(['a', 'b'], [['1', '2']])
    time.sleep(1.1)
    assert jobs._tsv_gz(['a', 'b'], [['1', '2']]) == first
