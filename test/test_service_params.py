"""Validation of /designs parameters.

These run without a bundle: `transcript_exists` is a stub, because the point is the
parsing rules rather than the reference data.
"""
import pytest
from starlette.datastructures import QueryParams

from service.params import ValidationError, parse_designs_query

KNOWN = 'ENST00000307102'


def parse(query, exists=lambda tid: True):
    return parse_designs_query(QueryParams(query), exists)


def test_bare_transcript_uses_the_engine_defaults():
    transcript, params = parse('transcript=' + KNOWN)
    assert transcript == KNOWN
    assert (params.pam, params.window, params.sg_len, params.edit) == \
        ('NGG', '4-8', 20, 'all')


def test_preset_resolves_through_the_engine_table():
    _, params = parse('transcript=%s&preset=ABE7.10' % KNOWN)
    assert (params.pam, params.window, params.edit) == ('NGG', '4-7', 'A-G')


def test_preset_still_takes_the_non_editor_overrides():
    _, params = parse('transcript=%s&preset=ABE7.10&intron_buffer=10&filter_gc=true'
                      % KNOWN)
    assert params.intron_buffer == 10 and params.filter_gc is True


def test_lowercase_input_is_normalised():
    transcript, params = parse('transcript=%s&pam=ngg' % KNOWN.lower())
    assert transcript == KNOWN and params.pam == 'NGG'


@pytest.mark.parametrize('query,expected', [
    ('', 'transcript is required'),
    ('transcript=../../etc/passwd', 'must look like'),
    ('transcript=ENST0000030710', 'must look like'),
    ('transcript=ENST123456789012', 'must look like'),
    ('transcript=%s&pam=XYZ' % KNOWN, 'IUPAC'),
    ('transcript=%s&pam=N' % KNOWN, 'IUPAC'),
    ('transcript=%s&preset=NOPE' % KNOWN, 'Unknown base editor'),
    ('transcript=%s&sg_len=16' % KNOWN, 'between 17 and 24'),
    ('transcript=%s&sg_len=25' % KNOWN, 'between 17 and 24'),
    ('transcript=%s&sg_len=x' % KNOWN, 'whole number'),
    ('transcript=%s&intron_buffer=101' % KNOWN, 'between 0 and 100'),
    ('transcript=%s&window=9-4' % KNOWN, 'start <= end'),
    ('transcript=%s&window=0-4' % KNOWN, 'start <= end'),
    ('transcript=%s&window=4-21' % KNOWN, 'start <= end'),
    ('transcript=%s&window=four' % KNOWN, 'look like 4-8'),
    ('transcript=%s&edit=DROP+TABLE' % KNOWN, 'edit must be one of'),
    ('transcript=%s&filter_gc=maybe' % KNOWN, 'true or false'),
])
def test_rejected(query, expected):
    with pytest.raises(ValidationError) as caught:
        parse(query)
    assert expected in str(caught.value)


def test_sql_injection_in_pam_is_a_validation_error_not_a_query():
    with pytest.raises(ValidationError):
        parse("transcript=%s&pam=' OR 1=1--" % KNOWN)


def test_unknown_parameters_are_refused_rather_than_ignored():
    # Ignoring them would let a user believe a setting applied to a result that
    # never saw it.
    with pytest.raises(ValidationError, match='Unknown parameter: evil'):
        parse('transcript=%s&evil=1' % KNOWN)


def test_a_repeated_parameter_is_refused():
    # Last-wins is how a request makes the validator and its reader disagree.
    with pytest.raises(ValidationError, match='more than once'):
        parse('transcript=%s&transcript=ENST00000000001' % KNOWN)


def test_preset_and_editor_parameters_together_are_refused():
    with pytest.raises(ValidationError, match='give one or the other'):
        parse('transcript=%s&preset=ABE7.10&pam=NGG' % KNOWN)


def test_an_overlong_value_is_refused_before_matching():
    with pytest.raises(ValidationError, match='too long'):
        parse('transcript=' + 'E' * 500)


def test_too_many_parameters_are_refused():
    query = '&'.join('p%d=1' % i for i in range(30))
    with pytest.raises(ValidationError, match='Too many parameters'):
        parse(query)


def test_a_transcript_missing_from_the_bundle_is_refused():
    with pytest.raises(ValidationError, match='not in the reference bundle'):
        parse('transcript=' + KNOWN, exists=lambda tid: False)


def test_the_bundle_is_not_consulted_for_a_malformed_id():
    """A junk ID must not reach the database at all."""
    asked = []
    with pytest.raises(ValidationError):
        parse('transcript=nonsense', exists=lambda tid: asked.append(tid) or True)
    assert asked == []
