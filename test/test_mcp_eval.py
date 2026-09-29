"""The evaluation: design tasks driven through the tools end to end.

Each task in mcp_eval_tasks.py is a question a researcher would ask, paired with the
route through the tool surface that should answer it. The runner replays that route,
checks every step was answerable, checks what the task says must be true of the
answers, and then checks the one thing that matters most:

    every guide handed back is a guide the engine actually designed, in its own text.

That last check is against test/golden, the same frozen output test_mcp_golden.py
compares byte for byte. Filtering, sorting and paging are ways of reading a result;
none of them may invent, alter or reorder a row's contents.

Every task needs a reference bundle, including the pasted-sequence one: the server
opens the bundle when it starts, even for a design that never consults it.
"""
import pytest

from mcp_eval_tasks import TASKS
from test_mcp_golden import read_golden

pytestmark = pytest.mark.bundle


def run_task(client, task):
    """Replays a task's route, returning each tool's answer by name."""
    results = {}
    for tool, build in task.plan:
        outcome = client.call_tool(tool, build(results))
        assert not outcome.is_error, '%s: %s failed: %s' % (
            task.name, tool, outcome.content)
        results[tool] = outcome.structured_content
    return results


@pytest.mark.parametrize('task', TASKS, ids=[t.name for t in TASKS])
def test_a_design_task_is_answerable_through_the_tools(task, client):
    """The route runs, says what the task expects, and invents nothing."""
    results = run_task(client, task)

    for label, check in task.expect:
        assert check(results), '%s: %s' % (task.name, label)

    page = results.get('query_guides')
    if page is None:
        return

    golden = read_golden(task.golden, 'sgrna_designs').splitlines()
    assert page['columns'] == golden[0].split('\t'), (
        '%s: the columns are not the file\'s columns' % task.name)
    invented = {'\t'.join(row) for row in page['rows']} - set(golden[1:])
    assert not invented, (
        '%s returned %d rows that are not in the frozen output: %s'
        % (task.name, len(invented), sorted(invented)[:2]))


# ---- the turns a model gets wrong first --------------------------------------

def test_an_unknown_editor_is_refused_rather_than_guessed(client):
    message = client.error('design_guides', {'transcript_id': 'ENST00000307102',
                                             'preset': 'BE9000'})
    assert 'BE9000' in message, 'the refusal should name what was not understood'


def test_an_ambiguous_gene_prefix_offers_the_choices(client):
    """The cue to ask the user, rather than to pick one."""
    resolved = client.call('resolve_gene', {'symbol': 'MAP2K'})
    assert resolved['gene'] == ''
    assert len(resolved['matches']) > 1
    assert 'MAP2K1' in resolved['matches']


def test_a_stale_handle_says_how_to_recover(client):
    assert 'design_guides again' in client.error('query_guides', {'run_id': 'a' * 64})
