"""The coverage panels: geometry, binning, the matrix, and the two filters.

Built from small hand-written frames rather than a design run, so a failure points at
the arithmetic instead of at the engine. The end-to-end checks against a real
transcript live in test_service_app.py.
"""
import numpy as np
import pandas as pd
import pytest

from service import coverage
from service.jobs import _tsv_gz
from service.params import TableView, ValidationError, parse_view_query
from service.results import ResultTable, UnknownFilterValue

# A two-exon gene on the plus strand, coding from 150.
PLUS = coverage.Geometry([(100, 200), (300, 400)], [(150, 200), (300, 350)], 1,
                         buffer=30)


# The columns the panels read.
COLUMNS = ['sgrna genomic position', 'sgRNA Strand', 'Mutation category',
           'Amino acid edits', 'Clinical significance']


def frame(rows):
    """A designs frame with only the columns the panels read."""
    return pd.DataFrame(rows, columns=COLUMNS)


def layout(geometry, *positions):
    return geometry.place(np.array(positions))[0].tolist()


def exon_of(geometry, *positions):
    return geometry.place(np.array(positions))[1].tolist()


def test_layout_runs_five_prime_to_three_prime():
    """Position along the drawing only ever increases with the coordinate."""
    places = layout(PLUS, *range(90, 410, 5))
    assert places == sorted(places)


def test_exons_are_drawn_to_scale_and_introns_are_not():
    """A 101 bp exon keeps its 101 units; the 99 bp intron between them does not."""
    start, end, next_start = layout(PLUS, 100, 200, 300)
    assert end - start == pytest.approx(100)
    intron = next_start - end
    assert intron < 99


def test_exon_parts_split_coding_from_untranslated():
    assert PLUS.parts(0) == [(100, 149, False), (150, 200, True)]
    assert PLUS.parts(1) == [(300, 350, True), (351, 400, False)]


def test_a_minus_strand_gene_is_laid_out_from_its_own_start():
    """Oriented coordinates put the transcript's first exon first.

    On the minus strand transcription begins at the high genomic coordinate, so the
    exon at 300-400 has to be drawn at the 5' end, not the 3'.
    """
    minus = coverage.Geometry([(100, 200), (300, 400)], [(150, 200)], -1, buffer=30)
    assert exon_of(minus, -400, -100) == [0, 1]
    first, last = layout(minus, -400, -100)
    assert first < last


def test_a_guide_is_drawn_as_the_worst_thing_it_does():
    """Two edits, one silent and one splice: the guide reads as the knockout."""
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Silent;Splice-donor', 'Lys5Lys;Exon1:+1', 'None;None'],
    ]), PLUS)
    assert table.class_counts['lof'] == 1
    assert table.class_counts.get('sil', 0) == 0


def test_a_guide_making_one_substitution_twice_counts_once():
    """The square says how many guides, because that is what clicking it returns."""
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense;Missense', 'Glu27Gly;Glu31Gly', 'None;None'],
    ]), PLUS)
    assert table.matrix_counts['Glu-Gly'] == 1
    assert table.matrix_sites['Glu-Gly'] == 2
    assert table.substitution_mask('Glu-Gly').sum() == 1


def test_splice_edits_stay_out_of_the_matrix():
    """'Exon2:+1' names no residue, so it fills no square."""
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Splice-donor', 'Exon1:+1', 'None'],
    ]), PLUS)
    assert not table.matrix_counts


def test_the_matrix_counts_pathogenic_recreations():
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense', 'Glu27Gly', 'Pathogenic'],
        ['165', 'sense', 'Missense', 'Glu27Gly', 'None'],
    ]), PLUS)
    assert table.matrix_counts['Glu-Gly'] == 2
    assert table.matrix_pathogenic['Glu-Gly'] == 1


def test_bins_are_bounded_whatever_the_guide_count():
    """A thousand guides draw no more bins than a dozen do."""
    rows = [[str(100 + i % 300), 'sense', 'Missense', 'Glu27Gly', 'None']
            for i in range(1000)]
    table = coverage.Coverage(frame(rows), PLUS)
    view = table.map_view()
    assert len(view['bins']) <= coverage.BINS
    assert len(view['exons']) == 2


def test_the_matrix_is_always_the_same_size():
    """441 squares whether or not the result fills any of them."""
    empty = coverage.Coverage(frame([]), PLUS).matrix_view()
    assert len(empty['cells']) == 21 * 21
    assert not empty['shaded']


def test_matrix_shading_spans_the_counts():
    """The fewest guides get the ramp's first colour and the most its last."""
    rows = ([['160', 'sense', 'Missense', 'Glu27Gly', 'None']] * 3
            + [['320', 'sense', 'Missense', 'Lys40Arg', 'None']])
    view = coverage.Coverage(frame(rows), PLUS).matrix_view()
    fills = {cell['key']: cell['fill'] for cell in view['cells'] if cell['fill']}
    assert fills == {'Glu-Gly': coverage.RAMP[-1], 'Lys-Arg': coverage.RAMP[0]}
    assert (view['low'], view['high']) == (1, 3)


def test_exon_mask_picks_out_one_exon():
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense', 'Glu27Gly', 'None'],
        ['320', 'sense', 'Missense', 'Lys40Arg', 'None'],
    ]), PLUS)
    assert table.exon_mask(1).tolist() == [True, False]
    assert table.exon_mask(2).tolist() == [False, True]
    assert table.has_exon(2)
    assert not table.has_exon(3)


def test_substitution_parsing_is_case_insensitive_and_bounded():
    assert coverage.parse_substitution('Glu-Gly') == ('Glu', 'Gly')
    assert coverage.parse_substitution('Xyz-Gly') is None
    assert coverage.parse_substitution('Glu-Gly-Ala') is None
    assert coverage.parse_substitution('') is None


def test_a_substitution_that_is_not_a_residue_pair_is_refused():
    with pytest.raises(ValidationError):
        parse_view_query(_query({'sub': 'drop table'}))
    with pytest.raises(ValidationError):
        parse_view_query(_query({'exon': '0'}))


def test_exon_and_substitution_survive_the_query_string():
    view = parse_view_query(_query({'exon': '3', 'sub': 'glu-gly'}))
    assert view.exon == 3
    assert view.sub == 'Glu-Gly'
    assert view.filtered


def _query(values):
    from starlette.datastructures import QueryParams
    return QueryParams(values)


def _designs_table(rows, geometry=PLUS):
    """A ResultTable over a frame, with its coverage attached the way the app does."""
    table = ResultTable(_tsv_gz(COLUMNS, rows), COLUMNS)
    table.coverage = coverage.Coverage(table.frame, geometry)
    return table


def test_the_table_answers_to_the_chart():
    table = _designs_table([
        ['160', 'sense', 'Missense', 'Glu27Gly', 'None'],
        ['165', 'sense', 'Silent', 'Lys5Lys', 'None'],
        ['320', 'sense', 'Missense', 'Glu40Gly', 'None'],
    ])
    assert table.select(TableView(sub='Glu-Gly')).matched == 2
    assert table.select(TableView(exon=1)).matched == 2
    # The two filters compose: exon 1 and Glu-Gly is the one guide in both.
    assert table.select(TableView(exon=1, sub='Glu-Gly')).matched == 1


def test_a_substitution_this_result_cannot_make_is_reported():
    """Not an empty table: without the matrix there was nothing to have clicked."""
    table = _designs_table([['160', 'sense', 'Missense', 'Glu27Gly', 'None']])
    with pytest.raises(UnknownFilterValue):
        table.select(TableView(sub='Trp-Ala'))
    with pytest.raises(UnknownFilterValue):
        table.select(TableView(exon=9))


def test_a_result_without_geometry_refuses_both_filters():
    """A transcript the bundle cannot place draws no panels, so neither filter is
    something a user could have clicked."""
    table = _designs_table([['160', 'sense', 'Missense', 'Glu27Gly', 'None']])
    table.coverage = None
    with pytest.raises(UnknownFilterValue):
        table.select(TableView(sub='Glu-Gly'))
    assert table.select(TableView()).matched == 1
