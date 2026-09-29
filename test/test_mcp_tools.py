"""How the tools behave at their edges.

The golden gate proves the tools reproduce the engine's output. This covers what the
gate cannot: how a bad argument comes back, what a handle does once its result is
gone, that a filtered page is still the engine's own rows, and that nothing the
engine prints can reach the protocol stream.
"""
import ast
import base64
import gzip
import json

import pytest

from bedesign import DESIGN_COLUMNS, DesignParams

from service.cachekey import manifest_key
from service.storage import LocalStorage

pytestmark = pytest.mark.bundle

MAP2K1 = 'ENST00000307102'
DESIGN = {'transcript_id': MAP2K1, 'pam': 'NGG', 'window': '4-8', 'sg_len': 20,
          'edit': 'all', 'intron_buffer': 30, 'filter_gc': False}


@pytest.fixture
def run(client):
    """One MAP2K1 design, shared by the tests that only read it back."""
    return client.call('design_guides', DESIGN)


# ---- the surface itself ------------------------------------------------------

def test_every_tool_is_advertised_as_read_only(client):
    listed = client.list_tools()
    names = {tool.name for tool in listed.tools}
    assert names == {'resolve_gene', 'list_transcripts', 'list_editors',
                     'design_guides', 'query_guides', 'get_clinvar_annotations',
                     'export_run'}
    for tool in listed.tools:
        assert tool.annotations.read_only_hint, tool.name
        assert tool.description, '%s has no description for the model' % tool.name


def test_list_editors_offers_only_editors_the_validator_accepts(client):
    """The dropdown and the allowlist are the same table, so this cannot drift."""
    editors = client.call('list_editors')
    assert editors['default'] == 'ABE8e-SpRY'
    for editor in editors['editors']:
        resolved = DesignParams.from_preset(editor['name'])
        assert editor['pam'] == resolved.pam
        assert editor['edit'] == resolved.edit
    assert editors['edits'] == ['C-T', 'A-G', 'all']


# ---- designing ---------------------------------------------------------------

def test_a_design_is_computed_once_and_addressed_by_its_parameters(client):
    """Asking the same question twice in a conversation costs one design."""
    first = client.call('design_guides', DESIGN)
    second = client.call('design_guides', DESIGN)
    assert first['run_id'] == second['run_id']
    # Different parameters are a different result, never a cache hit.
    other = client.call('design_guides', dict(DESIGN, edit='C-T'))
    assert other['run_id'] != first['run_id']


def test_a_design_reports_the_reference_versions_it_used(run):
    """A guide list is only meaningful against a stated reference."""
    assert 'Ensembl release' in run['reference']
    assert 'ClinVar' in run['clinvar']
    assert run['engine_version']


def test_giving_both_a_transcript_and_a_sequence_is_refused(client):
    message = client.error('design_guides',
                         {'transcript_id': MAP2K1, 'sequence': 'ACGT',
                          'sequence_name': 'x'})
    assert 'not both' in message


def test_giving_neither_a_transcript_nor_a_sequence_is_refused(client):
    assert 'either' in client.error('design_guides', {})


def test_a_sequence_needs_a_name(client):
    assert 'sequence_name' in client.error('design_guides',
                                         {'sequence': 'ACGTACGTACGTACGTACGTACGT'})


def test_a_sequence_that_is_not_nucleotides_is_refused(client):
    message = client.error('design_guides',
                         {'sequence': 'ACGTNNNNXYZ', 'sequence_name': 'junk'})
    assert 'nucleotides' in message


def test_an_unknown_transcript_says_so(client):
    message = client.error('design_guides', {'transcript_id': 'ENST99999999999'})
    assert 'not in the reference bundle' in message


def test_a_malformed_transcript_never_reaches_the_database(client):
    assert 'ENST' in client.error('design_guides', {'transcript_id': 'MAP2K1'})


def test_an_unknown_editor_is_refused_with_its_name(client):
    assert 'BE9000' in client.error('design_guides',
                                  {'transcript_id': MAP2K1, 'preset': 'BE9000'})


def test_a_preset_and_a_pam_together_are_refused(client):
    message = client.error('design_guides',
                         {'transcript_id': MAP2K1, 'preset': 'BE4max', 'pam': 'NGG'})
    assert 'preset already sets' in message


def test_an_intron_buffer_is_allowed_alongside_a_preset(client):
    """It is not part of the editor, so it is an override rather than a conflict."""
    designed = client.call('design_guides', {'transcript_id': MAP2K1,
                                             'preset': 'BE4max', 'intron_buffer': 0})
    assert designed['params']['intron_buffer'] == 0
    assert designed['params']['pam'] == DesignParams.from_preset('BE4max').pam


# ---- reading a run back ------------------------------------------------------

def test_a_page_is_fifty_rows_of_the_engines_own_text(run, client):
    page = client.call('query_guides', {'run_id': run['run_id']})
    assert page['columns'] == DESIGN_COLUMNS
    assert len(page['rows']) == 50
    assert page['matched'] == run['counts']['designs']
    assert all(isinstance(cell, str) for cell in page['rows'][0])


def test_a_page_past_the_end_lands_on_the_last_one(run, client):
    """A model that guessed a page number gets rows rather than an empty answer."""
    page = client.call('query_guides', {'run_id': run['run_id'], 'page': 9999})
    assert page['page'] == page['pages']
    assert page['rows']


def test_a_filter_narrows_to_a_subset_of_the_same_rows(run, client):
    everything = client.call('query_guides', {'run_id': run['run_id']})
    filtered = client.call('query_guides', {'run_id': run['run_id'],
                                            'consequence': 'lof'})
    assert 0 < filtered['matched'] < everything['matched']
    assert filtered['filters'] == {'consequence': 'lof'}


def test_columns_narrows_the_row_without_reordering_it(run, client):
    wanted = ['sgRNA sequence', 'Mutation category']
    page = client.call('query_guides', {'run_id': run['run_id'], 'columns': wanted})
    full = client.call('query_guides', {'run_id': run['run_id']})
    assert page['columns'] == wanted
    keep = [DESIGN_COLUMNS.index(name) for name in wanted]
    assert page['rows'][0] == [full['rows'][0][i] for i in keep]


def test_an_unknown_column_is_refused(run, client):
    assert 'Not a column' in client.error('query_guides',
                                        {'run_id': run['run_id'],
                                         'columns': ['sgRNA sequence', 'Nope']})


def test_a_significance_this_result_does_not_contain_is_refused(run, client):
    """Refused, not silently empty: the vocabulary is open, so an empty page would be
    indistinguishable from a typo."""
    message = client.error('query_guides',
                         {'run_id': run['run_id'], 'significance': 'Extremely Bad'})
    assert message


def test_sorting_reorders_the_same_rows(run, client):
    ascending = client.call('query_guides', {'run_id': run['run_id'],
                                             'sort': 'sgrna genomic position'})
    descending = client.call('query_guides', {'run_id': run['run_id'],
                                              'sort': 'sgrna genomic position',
                                              'dir': 'desc'})
    assert ascending['matched'] == descending['matched']
    assert ascending['rows'][0] != descending['rows'][0]


def test_a_handle_that_names_no_stored_result_says_what_to_do(client):
    message = client.error('query_guides', {'run_id': '0' * 64})
    assert 'design_guides again' in message


def test_a_handle_that_is_not_a_handle_is_refused(client):
    assert 'run_id' in client.error('query_guides', {'run_id': 'the last one'})


# ---- annotations and export --------------------------------------------------

def test_annotations_can_be_narrowed_to_one_guide(run, client):
    page = client.call('get_clinvar_annotations', {'run_id': run['run_id']})
    guide = page['rows'][0][0]
    only = client.call('get_clinvar_annotations',
                       {'run_id': run['run_id'], 'sgrna': guide})
    assert only['matched'] >= 1
    assert {row[0] for row in only['rows']} == {guide}


def test_matched_only_keeps_the_rows_that_found_a_variant(run, client):
    matched = client.call('get_clinvar_annotations',
                          {'run_id': run['run_id'], 'matched_only': True})
    assert matched['rows']
    assert all(len(row) == 28 for row in matched['rows'])


def test_export_names_the_file_the_cli_would_have_written(run, client):
    exported = client.call('export_run', {'run_id': run['run_id'], 'file': 'designs'})
    assert exported['filename'] == 'sgrna_designs_%s.txt' % MAP2K1
    assert exported['lines'] == run['counts']['designs'] + 1


def test_export_refuses_a_file_this_run_never_wrote(client):
    """A sequence run has no gene, so it has no annotations to export."""
    designed = client.call('design_guides',
                           {'sequence': 'ACGT' * 40, 'sequence_name': 'toy'})
    assert 'no gene' in client.error('export_run',
                                   {'run_id': designed['run_id'], 'file': 'clinvar'})


def test_an_oversize_text_export_names_the_other_encoding(run, client):
    """Never truncated: half a TSV that looks whole is worse than an error."""
    message = client.error('export_run',
                         {'run_id': run['run_id'], 'file': 'designs',
                          'max_bytes': 100})
    assert 'base64-gzip' in message


def test_the_gzip_encoding_carries_the_same_bytes(run, client):
    text = client.call('export_run', {'run_id': run['run_id'], 'file': 'designs'})
    packed = client.call('export_run', {'run_id': run['run_id'], 'file': 'designs',
                                        'encoding': 'base64-gzip'})
    unpacked = gzip.decompress(base64.b64decode(packed['data'])).decode('utf-8')
    # The stored file keeps the CRLF the CLI writes; the text export is what reading
    # that file back with universal newlines gives, which is what the goldens hold.
    assert unpacked.replace('\r\n', '\n') == text['data']
    assert packed['sha256'] == text['sha256']


# ---- sharing a store with the web service ----------------------------------

def test_a_run_the_web_service_designed_can_be_read_back(client, mcp_settings, run):
    """Both front doors write into one store under one key.

    service/jobs.py computes the same cache_key for the same transcript and
    parameters, and its manifest is the older, smaller shape. Reading it must serve
    the result that is sitting on disk rather than failing on a field the web service
    never wrote.
    """
    storage = LocalStorage(mcp_settings.results_dir)
    stored = json.loads(storage.get(manifest_key(run['run_id'])))
    # Exactly the fields service/jobs.py writes -- nothing this front door adds.
    storage.put(manifest_key(run['run_id']), json.dumps({
        'transcript_id': stored['transcript_id'], 'params': stored['params'],
        'engine_version': stored['engine_version'],
        'reference': stored['reference'], 'clinvar': stored['clinvar'],
        'designs': stored['designs'], 'errors': stored['errors'],
        'annotations': stored['annotations'],
        'runtime_seconds': stored['runtime_seconds'],
    }).encode())

    page = client.call('query_guides', {'run_id': run['run_id']})
    assert page['matched'] == stored['designs']
    exported = client.call('export_run', {'run_id': run['run_id'], 'file': 'designs'})
    assert exported['filename'] == 'sgrna_designs_%s.txt' % MAP2K1
    assert exported['sha256']
    annotations = client.call('get_clinvar_annotations', {'run_id': run['run_id']})
    assert annotations['rows']


# ---- the stdout hazard -------------------------------------------------------

def test_the_engine_never_writes_to_stdout():
    """On a stdio server, stdout is the JSON-RPC wire.

    A redirect cannot protect it: sys.stdout is process-global and the SDK answers
    several calls at once, so one call restoring it would uncover another's output.
    The engine and the transcript source therefore log their diagnostics instead of
    printing them, and this is the guard on that -- a print() added to either module
    would corrupt the protocol stream in a way no functional test would notice.
    """
    offenders = []
    for module in ('bedesign/engine.py', 'bedesign/transcript_source.py',
                   'mcp_server/runs.py', 'mcp_server/tools.py',
                   'mcp_server/server.py'):
        tree = ast.parse(open(module).read())
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == 'print'):
                offenders.append('%s:%d' % (module, node.lineno))
    assert not offenders, 'print() reaches the protocol stream at %s' % offenders


def test_a_design_writes_nothing_to_stdout(client, capsys):
    """The same invariant, exercised rather than read."""
    client.call('design_guides', DESIGN)
    assert capsys.readouterr().out == ''
