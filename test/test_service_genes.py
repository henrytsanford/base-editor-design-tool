"""Gene symbol lookup against a real bundle.

Marked `bundle` for the same reason the other reference tests are: the queries only
mean anything against the 646,577-row transcript table they were written for.
"""
import pytest

from service.references import References

pytestmark = pytest.mark.bundle

MAP2K1 = 'ENST00000307102'


@pytest.fixture(scope='module')
def references(refdata):
    return References(refdata)


def test_a_prefix_finds_the_gene_family(references):
    matches = references.search_genes('MAP2K')
    assert 'MAP2K1' in matches and 'MAP2K7' in matches
    assert all(name.startswith('MAP2K') for name in matches)


def test_the_prefix_search_is_index_backed(references, refdata):
    """Design doc 5.2: LIKE cannot use a BINARY-collation index, a range scan can.

    Asserted on the query plan rather than on a stopwatch, because a timing threshold
    would be flaky on a loaded machine while the plan is the thing that matters.
    """
    plan = references._db().execute(
        'EXPLAIN QUERY PLAN SELECT DISTINCT gene_name FROM transcript '
        'WHERE gene_name >= ? AND gene_name < ? ORDER BY gene_name LIMIT ?',
        ('MAP2K', 'MAP2K￿', 25)).fetchall()
    detail = ' '.join(row[-1] for row in plan)
    assert 'idx_transcript_gene_name' in detail
    assert 'SCAN' not in detail, 'the gene search must not scan the table'


def test_a_search_that_matches_nothing_is_empty_not_an_error(references):
    assert references.search_genes('ZZZZZZZZ') == []
    assert references.search_genes('') == []


def test_the_result_count_is_capped(references):
    """A one-letter query must not render every gene starting with it."""
    assert len(references.search_genes('A')) <= 25


def test_mane_select_leads_the_transcript_list(references):
    """Design doc 5: MANE Select is the preselected transcript."""
    transcripts = references.transcripts_for_gene('MAP2K1')
    assert transcripts[0]['transcript_id'] == MAP2K1
    assert transcripts[0]['mane_select'] == 1
    assert len(transcripts) > 1, 'MAP2K1 has several transcripts to choose between'


def test_a_transcript_row_carries_what_the_picker_shows(references):
    first = references.transcripts_for_gene('MAP2K1')[0]
    assert first['display_name'] == 'MAP2K1-201'
    assert first['biotype'] == 'protein_coding'
    assert first['seq_region'] == '15' and first['strand'] == 1
    assert first['start'] < first['end'] and first['cds_length'] > 0


def test_the_projection_leaves_the_sequences_behind(references):
    """cds_sequence and protein_sequence are why this database is 1.07 GB."""
    first = references.transcripts_for_gene('MAP2K1')[0]
    assert 'cds_sequence' not in first and 'protein_sequence' not in first


def test_an_unknown_gene_has_no_transcripts(references):
    assert references.transcripts_for_gene('NOTAGENE') == []
