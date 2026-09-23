"""LocalStorage, and the path safety it is responsible for."""
import os
import time

import pytest

from service.storage import BadKey, LocalStorage


@pytest.fixture
def store(tmp_path):
    return LocalStorage(str(tmp_path / 'results'))


def test_round_trip(store):
    store.put('results/ab/abcdef/designs.tsv.gz', b'payload')
    assert store.exists('results/ab/abcdef/designs.tsv.gz')
    assert store.get('results/ab/abcdef/designs.tsv.gz') == b'payload'


def test_a_missing_object_raises_key_error(store):
    assert not store.exists('results/ab/nothing/designs.tsv.gz')
    with pytest.raises(KeyError):
        store.get('results/ab/nothing/designs.tsv.gz')


def test_put_replaces_an_existing_object(store):
    store.put('results/ab/abcdef/manifest.json', b'old')
    store.put('results/ab/abcdef/manifest.json', b'new')
    assert store.get('results/ab/abcdef/manifest.json') == b'new'


def test_put_leaves_no_temporary_file_behind(store):
    store.put('results/ab/abcdef/manifest.json', b'x')
    written = os.listdir(os.path.join(store.root, 'results', 'ab', 'abcdef'))
    assert written == ['manifest.json']


@pytest.mark.parametrize('key', [
    '../escape', 'results/../../escape', '/etc/passwd', 'results//double',
    '', '   ', 'results/x/../../../../etc/passwd', ' results/a', 'results/a ',
    'results/./a', 'results/-leading/x',
])
def test_unsafe_keys_are_refused(store, key):
    with pytest.raises(BadKey):
        store.path(key)


def test_a_symlink_out_of_the_root_is_refused(store):
    """The segment rules cannot see this one, so the resolved path is checked too."""
    os.symlink('/etc', os.path.join(store.root, 'evil'))
    with pytest.raises(BadKey):
        store.path('evil/passwd')


def test_the_root_is_created_if_it_is_missing(tmp_path):
    root = tmp_path / 'deep' / 'nested'
    assert LocalStorage(str(root)).exists('anything') is False
    assert root.is_dir()


def put_result(store, digest, size, age):
    """One result group: a payload of `size` bytes and a manifest `age` s old."""
    prefix = 'results/%s/%s' % (digest[:2], digest)
    store.put(prefix + '/designs.tsv.gz', b'x' * size)
    store.put(prefix + '/manifest.json', b'{}')
    when = time.time() - age
    os.utime(store.path(prefix + '/manifest.json'), (when, when))
    return prefix


def test_evict_deletes_the_least_recently_used_results_down_to_the_budget(store):
    old = put_result(store, 'aa' * 32, 1000, age=300)
    middle = put_result(store, 'bb' * 32, 1000, age=200)
    new = put_result(store, 'cc' * 32, 1000, age=100)
    assert store.evict('results', 'manifest.json', 2100) == 1
    assert not store.exists(old + '/manifest.json')
    assert not os.path.exists(store.path(old))
    assert store.exists(middle + '/manifest.json') and store.exists(new + '/manifest.json')


def test_touching_a_result_keeps_it(store):
    """Eviction follows use, not creation: a popular result survives."""
    old = put_result(store, 'aa' * 32, 1000, age=300)
    newer = put_result(store, 'bb' * 32, 1000, age=100)
    store.touch(old + '/manifest.json')
    store.evict('results', 'manifest.json', 1500)
    assert store.exists(old + '/manifest.json')
    assert not store.exists(newer + '/manifest.json')


def test_evict_leaves_a_result_still_being_written(store):
    """No manifest means a job is still writing: not counted, not deleted."""
    store.put('results/dd/' + 'dd' * 32 + '/designs.tsv.gz', b'x' * 5000)
    put_result(store, 'aa' * 32, 1000, age=100)
    assert store.evict('results', 'manifest.json', 10) == 1
    assert store.exists('results/dd/' + 'dd' * 32 + '/designs.tsv.gz')


def test_evict_under_budget_deletes_nothing(store):
    put_result(store, 'aa' * 32, 1000, age=100)
    assert store.evict('results', 'manifest.json', 10 ** 6) == 0


def test_evict_clears_a_deletion_a_crash_left_behind(store):
    stale = os.path.join(store.root, 'results', 'aa', 'a' * 64 + '.evicting-1-0')
    os.makedirs(stale)
    open(os.path.join(stale, 'designs.tsv.gz'), 'wb').close()
    store.evict('results', 'manifest.json', 10 ** 6)
    assert not os.path.exists(stale)


def test_touching_a_missing_object_is_not_an_error(store):
    store.touch('results/ab/gone/manifest.json')
