"""The job table: deduplication, the busy paths, failures, crashes and the timeout.

The pool is replaced with a fake so these run offline and deterministically. What is
under test is the bookkeeping in JobPool, not the design engine.
"""
import threading
import time
from concurrent.futures import Future
from concurrent.futures.process import BrokenProcessPool

import pytest

import service.jobs as jobs
from service.config import Settings
from service.jobs import ClientBusy, JobPool, JobTimeout, PoolBusy


class FakeExecutor(object):
    def __init__(self, **kwargs):
        self.calls = []
        self.broken = False
        self.shut_down = False
        # ProcessPoolExecutor fails a broken pool's futures while holding this
        # non-reentrant lock, which is what makes shutdown() from a done-callback
        # deadlock. Modelled so the tests below can hold it the same way.
        self.shutdown_lock = threading.Lock()

    def submit(self, fn, *args):
        if self.broken:
            raise BrokenProcessPool('a child process terminated abruptly')
        future = Future()
        self.calls.append((args, future))
        return future

    def shutdown(self, **kwargs):
        if not self.shutdown_lock.acquire(timeout=1):
            raise AssertionError('shutdown() would deadlock: its lock is held')
        self.shutdown_lock.release()
        self.shut_down = True

    def break_(self):
        """Fails every pending future the way a broken pool does: under its lock."""
        self.broken = True
        with self.shutdown_lock:
            for _, future in self.calls:
                if not future.done():
                    future.set_exception(BrokenProcessPool('terminated abruptly'))


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


def test_a_client_runs_one_job_at_a_time(pool):
    """Without this, one client could hold every worker inside its rate limit."""
    assert pool.submit('a', 'ENST1', None, client='1.2.3.4') is True
    assert pool.client_busy('1.2.3.4') is True
    with pytest.raises(ClientBusy):
        pool.submit('b', 'ENST2', None, client='1.2.3.4')
    assert pool.submit('c', 'ENST3', None, client='5.6.7.8') is True


def test_joining_a_running_job_is_not_a_second_job(pool):
    """A client polling its own design, or anyone asking for one already running,
    joins it rather than being turned away."""
    pool.submit('a', 'ENST1', None, client='1.2.3.4')
    assert pool.submit('a', 'ENST1', None, client='1.2.3.4') is False
    assert pool.submit('a', 'ENST1', None, client='5.6.7.8') is False
    assert pool.client_busy('5.6.7.8') is False


def test_a_client_is_free_once_its_job_ends(pool):
    pool.submit('a', 'ENST1', None, client='1.2.3.4')
    pool._pool.calls[0][1].set_exception(RuntimeError('boom'))
    assert pool.client_busy('1.2.3.4') is False
    assert pool.submit('b', 'ENST2', None, client='1.2.3.4') is True


def test_a_dead_worker_replaces_the_pool(pool):
    """ProcessPoolExecutor never recovers from a worker killed outright; every
    later submit would raise, and the service would be down until restarted."""
    pool.submit('a', 'ENST1', None, client='1.2.3.4')
    first = pool._pool
    first.break_()
    assert pool._pool is not first
    assert pool.failure('a') == 'This design failed. Reload to try again.'
    assert pool.client_busy('1.2.3.4') is False
    assert pool.submit('b', 'ENST2', None) is True
    assert len(pool._pool.calls) == 1


def test_jobs_failing_together_replace_the_pool_once(pool):
    pool.submit('a', 'ENST1', None)
    pool.submit('b', 'ENST2', None)
    first = pool._pool
    first.break_()
    second = pool._pool
    assert second is not first
    assert pool.submit('c', 'ENST3', None) is True
    assert pool._pool is second and len(second.calls) == 1


def test_replacing_a_broken_pool_does_not_shut_it_down_from_its_callback(pool):
    """The deadlock a real SIGKILL found: ProcessPoolExecutor runs done-callbacks
    while holding its shutdown lock, so shutdown() there never returns and every
    request then blocks on JobPool's lock."""
    pool.submit('a', 'ENST1', None)
    first = pool._pool
    first.break_()  # the fake's shutdown() raises if called while this lock is held
    assert not first.shut_down


def test_a_pool_that_broke_while_idle_is_replaced_on_submit(pool):
    """A worker can die with no job running, so no future reports it."""
    first = pool._pool
    first.broken = True
    assert pool.submit('a', 'ENST1', None) is True
    assert pool._pool is not first and len(pool._pool.calls) == 1


def test_health_reports_a_crash_until_a_job_next_finishes(pool):
    assert pool.healthy() is True
    pool.submit('a', 'ENST1', None)
    pool._pool.break_()
    assert pool.healthy() is False
    pool.submit('b', 'ENST2', None)
    pool._pool.calls[0][1].set_result({})
    assert pool.healthy() is True


def test_an_ordinary_failure_is_not_a_crash(pool):
    pool.submit('a', 'ENST1', None)
    first = pool._pool
    first.calls[0][1].set_exception(JobTimeout('too slow'))
    assert pool.healthy() is True and pool._pool is first


def hooked_pool(monkeypatch, hook):
    monkeypatch.setattr(jobs, 'ProcessPoolExecutor',
                        lambda **kwargs: FakeExecutor(**kwargs))
    return JobPool(FakeReferences(), Settings(max_jobs=2), after_job=hook)


def test_the_after_job_hook_runs_once_a_result_is_written(monkeypatch):
    """What keeps RESULTS_DIR under its budget: it runs after each new result."""
    calls = []
    pool = hooked_pool(monkeypatch, lambda: calls.append(1))
    pool.submit('a', 'ENST1', None)
    pool._pool.calls[0][1].set_result({})
    assert calls == [1]


def test_the_after_job_hook_does_not_run_for_a_failure(monkeypatch):
    calls = []
    pool = hooked_pool(monkeypatch, lambda: calls.append(1))
    pool.submit('a', 'ENST1', None)
    pool._pool.calls[0][1].set_exception(RuntimeError('boom'))
    assert calls == []


def test_a_failing_hook_does_not_fail_the_design(monkeypatch):
    def hook():
        raise OSError('disk trouble')

    pool = hooked_pool(monkeypatch, hook)
    pool.submit('a', 'ENST1', None, client='1.2.3.4')
    pool._pool.calls[0][1].set_result({})
    assert pool.failure('a') is None and pool.running('a') is False
    assert pool.client_busy('1.2.3.4') is False


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
