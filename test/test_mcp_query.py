"""The Query adapter, and the rule that MCP inherits the web app's validation.

The point of the adapter is that service/params.py validates tool arguments without
being changed, so what is checked here is mostly that the inheritance is real: a rule
written for the web app applies to the tools, and an argument the model left unset is
read as absent rather than as an empty string.

Offline: no bundle, no reference data, no MCP.
"""
import pytest

from bedesign.engine import DesignParams

from mcp_server.query import Query
from service import params
from service.params import ValidationError


def test_an_unset_argument_is_absent_rather_than_empty():
    """None means the model did not give the argument at all.

    This is what lets a preset coexist with the nine other editor arguments sitting
    at their defaults: params reads an absent value, not a conflicting one.
    """
    query = Query({'preset': 'BE4max', 'pam': None, 'window': None, 'sg_len': None})
    assert set(query.keys()) == {'preset'}
    assert query.getlist('pam') == []
    assert query.getlist('preset') == ['BE4max']


def test_booleans_arrive_as_the_words_params_reads():
    assert Query({'filter_gc': True}).getlist('filter_gc') == ['true']
    assert Query({'filter_gc': False}).getlist('filter_gc') == ['false']
    assert params.parse_editor_params(Query({'filter_gc': True})).filter_gc is True


def test_numbers_arrive_as_text():
    assert Query({'sg_len': 21}).getlist('sg_len') == ['21']
    assert params.parse_editor_params(Query({'sg_len': 21})).sg_len == 21


def test_a_preset_resolves_through_the_engines_own_table():
    resolved = params.parse_editor_params(Query({'preset': 'ABE7.10'}))
    assert (resolved.pam, resolved.window, resolved.edit) == ('NGG', '4-7', 'A-G')


def test_a_preset_alongside_an_editor_parameter_is_refused():
    """Which one would win? The web app refuses this, so the tools do too."""
    with pytest.raises(ValidationError) as caught:
        params.parse_editor_params(Query({'preset': 'BE4max', 'pam': 'NGG'}))
    assert 'preset already sets pam' in str(caught.value)


def test_an_unknown_editor_is_refused():
    with pytest.raises(ValidationError):
        params.parse_editor_params(Query({'preset': 'BE9000'}))


@pytest.mark.parametrize('bad', [
    {'pam': 'NXX'},            # X is not an IUPAC code
    {'pam': 'N'},              # shorter than the 2-8 the pattern allows
    {'sg_len': 30},            # past SG_LEN_RANGE
    {'sg_len': 3},
    {'intron_buffer': 500},    # past INTRON_BUFFER_RANGE
    {'window': '8-4'},         # end before start
    {'window': '4-40'},        # past sg_len
    {'window': 'four-eight'},
    {'edit': 'C-G'},           # not a deaminase this engine has
])
def test_the_web_apps_editor_rules_apply_to_tool_arguments(bad):
    with pytest.raises(ValidationError):
        params.parse_editor_params(Query(bad))


def test_the_window_is_checked_against_the_sg_len_the_request_uses():
    """Not against the default -- a longer guide legitimately allows a later window."""
    assert params.parse_editor_params(Query({'sg_len': 24, 'window': '4-24'})).window == '4-24'
    with pytest.raises(ValidationError):
        params.parse_editor_params(Query({'sg_len': 20, 'window': '4-24'}))


def test_defaults_come_from_the_engine():
    """An empty request is DesignParams(), so the tools cannot drift from the CLI."""
    assert params.parse_editor_params(Query({})) == DesignParams()


def test_a_view_is_parsed_from_named_arguments():
    view = params.parse_view_query(Query({
        'consequence': 'mis', 'significance': 'Pathogenic', 'page': 3,
        'hide_bsmbi': True, 'exon': 5, 'dir': 'desc', 'sort': 'sgRNA sequence'}))
    assert (view.consequence, view.significance, view.page) == ('mis', 'Pathogenic', 3)
    assert (view.hide_bsmbi, view.exon, view.dir) == (True, 5, 'desc')


def test_an_unset_filter_is_not_a_filter():
    view = params.parse_view_query(Query({'consequence': None, 'exon': None}))
    assert not view.filtered


@pytest.mark.parametrize('bad', [
    {'consequence': 'nonsense'},          # not one of the five classes
    {'strand': 'forward'},                # sense / antisense
    {'deaminase': 'all'},                 # a design choice, not a view one
    {'sort': 'Not A Column'},
    {'dir': 'sideways'},
    {'sub': 'Glu>Gly'},                   # the matrix key is Glu-Gly
    {'exon': 0},
])
def test_the_web_apps_view_rules_apply_too(bad):
    with pytest.raises(ValidationError):
        params.parse_view_query(Query(bad))
