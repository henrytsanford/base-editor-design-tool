"""LocalStorage, and the path safety it is responsible for."""
import os

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
