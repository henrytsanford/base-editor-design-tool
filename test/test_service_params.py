"""Validation of query parameters, for every route that takes them.

These run without a bundle: `transcript_exists` is a stub, because the point is the
parsing rules rather than the reference data.
"""
import pytest
from starlette.datastructures import QueryParams

from service.params import (DESIGN_PARAMS, VIEW_PARAMS, ValidationError,
                            parse_designs_query, parse_download_query,
                            parse_genes_query, parse_view_query)

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


# --- The view half: filters, sort and paging (design doc 5, slice 2) -------------


def view(query):
    return parse_view_query(QueryParams(query))


def test_a_view_parameter_can_never_shadow_a_design_one():
    """The two halves are disjoint, or a filter would change the cache key.

    'deaminase' rather than 'edit' is the whole reason this holds: both name the same
    C-T/A-G vocabulary, but only one of them describes what to compute.
    """
    assert not set(VIEW_PARAMS) & set(DESIGN_PARAMS)
    assert 'edit' in DESIGN_PARAMS and 'deaminase' in VIEW_PARAMS


def test_filters_do_not_change_the_cache_key():
    """Filtering re-renders a cached result; it must never recompute one."""
    plain = parse('transcript=%s&preset=ABE7.10' % KNOWN)
    filtered = parse('transcript=%s&preset=ABE7.10&mutation=Missense&sort=PAM&page=3'
                     % KNOWN)
    assert plain == filtered


def test_a_filtered_url_is_a_legal_designs_request():
    transcript, _ = parse('transcript=%s&mutation=Missense&hide_4t=true' % KNOWN)
    assert transcript == KNOWN


def test_view_defaults_show_the_unfiltered_first_page():
    parsed = view('')
    assert parsed.page == 1 and parsed.dir == 'asc' and parsed.sort == ''
    assert not parsed.filtered


def test_any_filter_marks_the_view_as_filtered():
    assert view('mutation=Missense').filtered
    assert view('hide_4t=true').filtered
    assert not view('sort=PAM&page=2').filtered, 'sorting is not filtering'


def test_sort_must_name_a_column_of_the_table():
    assert view('sort=%23 edits').sort == '# edits'
    with pytest.raises(ValidationError):
        view('sort=DROP TABLE')


@pytest.mark.parametrize('query', [
    'dir=sideways',
    'strand=both',
    'deaminase=G-C',
    'page=0',
    'page=nine',
    'mutation=<script>',
    'hide_4t=perhaps',
])
def test_a_bad_view_value_is_refused(query):
    with pytest.raises(ValidationError):
        view(query)


def test_a_repeated_view_parameter_is_refused_like_a_design_one():
    with pytest.raises(ValidationError):
        view('page=1&page=2')


def test_an_open_vocabulary_filter_is_checked_for_shape_only():
    """ClinVar's classification list grows without asking us, so the value is matched
    against the result's own vocabulary later -- see test_service_results."""
    assert view('significance=Benign/Likely benign').significance == \
        'Benign/Likely benign'
    assert view('significance=Conflicting classifications of pathogenicity')


def test_a_blank_field_means_the_user_left_it_alone():
    """A GET form submits every field it has, including the empty ones."""
    _, params = parse('transcript=%s&preset=ABE7.10&pam=&window=&sg_len=' % KNOWN)
    assert params.pam == 'NGG' and params.window == '4-7'
    assert view('mutation=&sort=&page=').mutation == ''


# --- /genes and /designs/download -----------------------------------------------


def test_a_gene_query_carries_the_editor_choice():
    symbol, params = parse_genes_query(QueryParams('q=map2k1&preset=ABE7.10'))
    assert symbol == 'MAP2K1' and params.edit == 'A-G'


def test_an_empty_search_box_is_not_an_error():
    symbol, _ = parse_genes_query(QueryParams('q='))
    assert symbol == ''


@pytest.mark.parametrize('query', ['q=<script>', 'q=' + 'A' * 40, 'q=MAP2K1&evil=1'])
def test_a_bad_gene_query_is_refused(query):
    with pytest.raises(ValidationError):
        parse_genes_query(QueryParams(query))


def test_a_download_names_a_file_from_the_fixed_table():
    _, _, name = parse_download_query(
        QueryParams('transcript=%s&file=clinvar' % KNOWN), lambda tid: True)
    assert name == 'clinvar'


def test_a_download_defaults_to_the_designs_file():
    _, _, name = parse_download_query(
        QueryParams('transcript=' + KNOWN), lambda tid: True)
    assert name == 'designs'


@pytest.mark.parametrize('query', [
    'transcript=%s&file=../../etc/passwd' % KNOWN,
    'transcript=%s&file=manifest' % KNOWN,
    # Filters are refused here: a download is always the complete result, so it can
    # never hand back a partial file that looks like the whole one.
    'transcript=%s&file=designs&mutation=Missense' % KNOWN,
])
def test_a_bad_download_request_is_refused(query):
    with pytest.raises(ValidationError):
        parse_download_query(QueryParams(query), lambda tid: True)
