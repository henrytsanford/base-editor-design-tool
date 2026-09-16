"""The cache key (design doc 3.1)."""
import json

import pytest

from bedesign import DesignParams
from service.cachekey import (RESULT_FILES, cache_key, key_fields, manifest_key,
                              result_key, result_prefix)

ARGS = ('ENST00000307102', DesignParams(), '116', '2026-09-09', '3.0.0')


def test_a_key_is_a_sha256_digest():
    key = cache_key(*ARGS)
    assert len(key) == 64 and all(c in '0123456789abcdef' for c in key)


def test_the_same_request_gives_the_same_key():
    assert cache_key(*ARGS) == cache_key(*ARGS)


def test_a_preset_and_the_parameters_it_resolves_to_share_a_result():
    """Design doc 3.1: identical resolved values, one cached result."""
    preset = DesignParams.from_preset('ABE7.10')
    explicit = DesignParams(pam='NGG', window='4-7', sg_len=20, edit='A-G')
    assert cache_key('ENST1', preset, '116', 'v', '3.0.0') == \
        cache_key('ENST1', explicit, '116', 'v', '3.0.0')


@pytest.mark.parametrize('field,value', [
    ('pam', 'NGA'), ('window', '4-7'), ('sg_len', 21), ('edit', 'C-T'),
    ('intron_buffer', 20), ('filter_gc', True),
])
def test_every_design_parameter_changes_the_key(field, value):
    other = DesignParams(**{field: value})
    assert cache_key('ENST1', other, '116', 'v', '3.0.0') != \
        cache_key('ENST1', DesignParams(), '116', 'v', '3.0.0')


@pytest.mark.parametrize('index,label', [(0, 'transcript'), (2, 'release'),
                                         (3, 'clinvar'), (4, 'engine_version')])
def test_every_context_field_changes_the_key(index, label):
    changed = list(ARGS)
    changed[index] = 'different'
    assert cache_key(*changed) != cache_key(*ARGS)


def test_the_key_covers_exactly_the_documented_fields():
    assert set(key_fields(*ARGS)) == {
        'transcript_id', 'ensembl_release', 'clinvar_version', 'pam', 'window',
        'sg_len', 'edit', 'intron_buffer', 'filter_gc', 'engine_version'}


def test_the_key_does_not_depend_on_field_order():
    fields = key_fields(*ARGS)
    forward = json.dumps(fields, sort_keys=True, separators=(',', ':'))
    backward = json.dumps(dict(reversed(list(fields.items()))), sort_keys=True,
                          separators=(',', ':'))
    assert forward == backward


def test_object_keys_contain_only_the_digest_and_fixed_names():
    """The digest is the only thing that becomes a path segment."""
    key = cache_key(*ARGS)
    assert result_prefix(key) == 'results/%s/%s' % (key[:2], key)
    for name, filename in RESULT_FILES.items():
        assert result_key(key, name).endswith('/' + filename)
    assert manifest_key(key).endswith('/manifest.json')


def test_an_unknown_result_name_raises_rather_than_building_a_path():
    with pytest.raises(KeyError):
        result_key(cache_key(*ARGS), '../../etc/passwd')
