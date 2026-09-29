"""The MCP surface, against the frozen engine output.

The risk in putting a conversational layer over a scientific engine is that the
layer quietly changes the science. This is what rules that out: the same six cases
test_golden.py runs through the CLI are run through the MCP tools, and every output
file is compared byte for byte against the same fixtures.

If a tool can reproduce every byte the CLI writes, the conversational surface is a
way of reading the engine's output rather than a second implementation of it.

The panel and the fixture reader are imported from test_golden rather than repeated,
so the two suites cannot come to disagree about what the golden output is.
"""
from Bio import SeqIO

import pytest

from bedesign import ANNOTATION_COLUMNS

# pytest puts test/ on sys.path, so the golden panel is importable by name.
from test_golden import NUC_OUTPUTS, PANEL, TID_OUTPUTS, read_golden

# What run_design passes the CLI for a golden run: its fixed flags, plus --edit all.
GOLDEN_PARAMS = {'pam': 'NGG', 'window': '4-8', 'sg_len': 20, 'edit': 'all',
                 'intron_buffer': 30, 'filter_gc': False}

# The fixture name each stored file answers to.
FILES = {'sgrna_designs': 'designs', 'error_report': 'errors',
         'clinvar_annotations': 'clinvar'}

GFP_FASTA = 'Sample_data/GFP.fasta'


def gfp_sequence():
    """The GFP sequence as the CLI reads it.

    Through SeqIO, the way read_args does, because case is load bearing: the
    nucleotide path marks intronic bases by lower case and GFP.fasta has three.
    """
    with open(GFP_FASTA) as fh:
        record = next(SeqIO.parse(fh, 'fasta'))
    return record.id, str(record.seq)


def check_exports(client, run_id, name, outputs):
    """Every file this run wrote, against its fixture."""
    for output, _template in outputs:
        exported = client.call('export_run', {'run_id': run_id,
                                              'file': FILES[output],
                                              'encoding': 'text'})
        assert exported['data'] == read_golden(name, output), (
            '%s/%s differs from its golden fixture. The MCP server and the CLI '
            'must write the same bytes.' % (name, output))


@pytest.mark.bundle
@pytest.mark.parametrize('name,transcript,gene', PANEL, ids=[p[0] for p in PANEL])
def test_the_tools_reproduce_the_golden_output(client, name, transcript, gene):
    """A design reached the way a conversation reaches it, byte for byte."""
    resolved = client.call('resolve_gene', {'symbol': gene})
    assert resolved['gene'] == gene

    listed = client.call('list_transcripts', {'gene': gene})
    assert transcript in [t['transcript_id'] for t in listed['transcripts']]

    run = client.call('design_guides', dict(GOLDEN_PARAMS, transcript_id=transcript))
    check_exports(client, run['run_id'], name, TID_OUTPUTS)

    # The summary counts what the file contains, so a count cannot drift from it.
    assert run['counts']['designs'] == len(
        read_golden(name, 'sgrna_designs').splitlines()) - 1
    assert run['counts']['annotations'] == len(
        read_golden(name, 'clinvar_annotations').splitlines()) - 1


@pytest.mark.bundle
def test_the_nucleotide_path_reproduces_its_golden_output(client):
    """A pasted sequence, which needs no gene and writes no annotations.

    Bundle-marked even so: the design consults no reference, but the server opens one
    when it starts, so the fixture cannot build without it.
    """
    name, sequence = gfp_sequence()
    run = client.call('design_guides',
                      dict(GOLDEN_PARAMS, sequence=sequence, sequence_name=name))
    assert run['files'] == ['designs', 'errors']
    check_exports(client, run['run_id'], 'gfp', NUC_OUTPUTS)


@pytest.mark.bundle
def test_the_first_page_of_rows_is_the_first_page_of_the_fixture(client):
    """query_guides pages the rows export_run writes, in the file's own order."""
    run = client.call('design_guides',
                      dict(GOLDEN_PARAMS, transcript_id='ENST00000307102'))
    page = client.call('query_guides', {'run_id': run['run_id'], 'page': 1})
    golden = read_golden('map2k1', 'sgrna_designs').splitlines()
    assert page['columns'] == golden[0].split('\t')
    assert ['\t'.join(row) for row in page['rows']] == golden[1:51]


@pytest.mark.bundle
def test_annotation_rows_keep_their_natural_length(client):
    """28 fields on a ClinVar match, 17 otherwise -- never padded to the header.

    A whole-file comparison would catch padding too, but only as a wall of diff.
    This says which invariant broke.
    """
    run = client.call('design_guides',
                      dict(GOLDEN_PARAMS, transcript_id='ENST00000307102'))
    annotations = client.call('get_clinvar_annotations', {'run_id': run['run_id']})
    lengths = {len(row) for row in annotations['rows']}
    assert lengths <= {17, len(ANNOTATION_COLUMNS)}, lengths
    assert 17 in lengths, 'every row matched a SNP; the ragged case is not covered'
    assert annotations['ragged'] is True


@pytest.mark.bundle
def test_a_non_coding_transcript_reports_an_error_rather_than_failing(client):
    """MALAT1 has no CDS. The engine writes an error row; the tool must say so."""
    run = client.call('design_guides',
                      dict(GOLDEN_PARAMS, transcript_id='ENST00000620902'))
    assert run['counts']['designs'] == 0
    assert run['counts']['errors'] == 1
    check_exports(client, run['run_id'], 'malat1', TID_OUTPUTS)
