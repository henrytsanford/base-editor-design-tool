"""The coverage panels: geometry, binning, the matrix, and the two filters.

Built from small hand-written frames rather than a design run, so a failure points at
the arithmetic instead of at the engine. The end-to-end checks against a real
transcript live in test_service_app.py.
"""
import numpy as np
import pytest

from service import coverage
from service.jobs import _tsv_gz
from service.params import TableView
from service.results import ResultTable, UnknownFilterValue

# A two-exon gene on the plus strand, coding from 150.
PLUS = coverage.Geometry([(100, 200), (300, 400)], [(150, 200), (300, 350)], 1,
                         buffer=30)


# The columns the panels read.
COLUMNS = ['sgrna genomic position', 'sgRNA Strand', 'Mutation category',
           'Amino acid edits', 'Clinical significance']


def frame(rows):
    """A designs table with only the columns the panels read."""
    return ResultTable(_tsv_gz(COLUMNS, rows), COLUMNS)


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


# Coding 1000-1100, between a 5 kb 5' UTR and a 3 kb 3' UTR on the same exon, then a
# second exon that is all UTR.
LONG_UTR = coverage.Geometry([(-4000, 4100), (9000, 9500)], [(1000, 1100)], 1,
                             buffer=30)


def test_a_long_utr_keeps_only_the_buffer_beside_the_coding_sequence():
    """The 5 kb UTR draws as the buffer either side of a stub, like an intron;
    the coding sequence and the buffer next to it still run base for base."""
    utr_start, cds_start, cds_end = layout(LONG_UTR, -4000, 1000, 1101)
    assert cds_end - cds_start == pytest.approx(101)
    assert cds_start - utr_start == pytest.approx(2 * 30 + coverage.STUB)
    near, edge = layout(LONG_UTR, 970, 1000)
    assert edge - near == pytest.approx(30)


def test_a_utr_only_exon_is_squeezed_but_still_drawn():
    start, end = layout(LONG_UTR, 9000, 9501)
    assert end - start == pytest.approx(2 * 30 + coverage.STUB)
    assert LONG_UTR.width < 700


def test_an_exon_number_sits_on_its_coding_part():
    """Exon 1 is UTR from 100 and coding from 150: the number centres on 150-200,
    not on the whole exon, where it would half sit on the thin UTR block."""
    exon = coverage.Coverage(frame([]), PLUS).map_view()['exons'][0]
    scale = (coverage.WIDTH - 2 * coverage.MARGIN) / PLUS.width
    start, end = layout(PLUS, 150, 201)
    assert exon['label_x'] == pytest.approx(
        coverage.MARGIN + (start + end) / 2 * scale, abs=0.01)


def test_exon_parts_split_coding_from_untranslated():
    assert PLUS.parts[0] == [(100, 149, False), (150, 200, True)]
    assert PLUS.parts[1] == [(300, 350, True), (351, 400, False)]


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


def test_a_selected_substitution_fades_the_bins_without_it():
    """Exon 1's bin holds the Glu-Gly guide and keeps its colour; exon 2's does not."""
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense', 'Glu27Gly', 'None'],
        ['380', 'sense', 'Missense', 'Lys40Arg', 'None'],
    ]), PLUS)
    assert [b['dim'] for b in table.map_view(selected_substitution='Glu-Gly')['bins']] \
        == [False, True]
    assert not any(b['dim'] for b in table.map_view()['bins'])


def test_a_selected_substitution_lights_only_its_own_guides():
    """In a bin shared with other guides only the selected ones stay lit, drawn in
    the square's colour: this guide's worst edit is missense, but it was picked
    from the silent square, so it lights green."""
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense;Silent', 'Glu27Gly;Leu28Leu', 'None;None'],
        ['160', 'sense', 'Missense', 'Lys40Arg', 'None'],
    ]), PLUS)
    [bin_] = table.map_view(selected_substitution='Leu-Leu')['bins']
    lit = [(seg['cls'], seg['height']) for seg in bin_['segments'] if not seg['dim']]
    faded = [seg['cls'] for seg in bin_['segments'] if seg['dim']]
    assert [cls for cls, _ in lit] == ['sil'] and faded == ['mis']
    # the two halves still stack to the bin's full height
    assert sum(seg['height'] for seg in bin_['segments']) == coverage.BAR_MAX
    assert not bin_['dim']


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
    """The fewest guides get the palest fill and the most a solid one."""
    rows = ([['160', 'sense', 'Missense', 'Glu27Gly', 'None']] * 3
            + [['320', 'sense', 'Missense', 'Lys40Arg', 'None']])
    view = coverage.Coverage(frame(rows), PLUS).matrix_view()
    shades = {cell['key']: cell['opacity'] for cell in view['cells'] if cell['count']}
    assert shades == {'Glu-Gly': 1.0, 'Lys-Arg': coverage.MIN_OPACITY}
    assert (view['low'], view['high']) == (1, 3)
    assert view['steps'][0] == coverage.MIN_OPACITY and view['steps'][-1] == 1.0


def test_matrix_squares_take_the_map_consequence_colours():
    """A stop is loss of function and the diagonal silent, as on the map, and the
    legend lists only the classes the result fills."""
    rows = [['160', 'sense', 'Missense', 'Glu27Gly', 'None'],
            ['170', 'sense', 'Nonsense', 'Gln28Ter', 'None'],
            ['180', 'sense', 'Silent', 'Leu29Leu', 'None']]
    view = coverage.Coverage(frame(rows), PLUS).matrix_view()
    classes = {cell['key']: cell['cls'] for cell in view['cells'] if cell['count']}
    assert classes == {'Glu-Gly': 'mis', 'Gln-Ter': 'lof', 'Leu-Leu': 'sil'}
    assert [entry['cls'] for entry in view['classes']] == ['lof', 'mis', 'sil']
    assert all(cell['cls'] == 'off' for cell in view['cells'] if not cell['count'])


def test_exon_mask_picks_out_one_exon():
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense', 'Glu27Gly', 'None'],
        ['320', 'sense', 'Missense', 'Lys40Arg', 'None'],
    ]), PLUS)
    assert table.exon_mask(1).tolist() == [True, False]
    assert table.exon_mask(2).tolist() == [False, True]
    assert table.has_exon(2)
    assert not table.has_exon(3)


def test_the_form_offers_only_exons_with_guides():
    """Unless a link already selected the empty one: the form must still show it."""
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense', 'Glu27Gly', 'None'],
        ['170', 'sense', 'Missense', 'Lys40Arg', 'None'],
    ]), PLUS)
    assert table.exon_options() == [(1, 2)]
    assert table.exon_options(selected=2) == [(1, 2), (2, 0)]


def test_the_form_offers_the_squares_the_matrix_links_in_its_order():
    table = coverage.Coverage(frame([
        ['160', 'sense', 'Missense', 'Glu27Gly', 'None'],
        ['170', 'sense', 'Missense;Nonsense', 'Lys40Arg;Glu41Ter', 'None;None'],
        ['180', 'sense', 'Missense', 'Lys42Arg', 'None'],
    ]), PLUS)
    options = table.substitution_options()
    # Lys before Glu: the matrix groups positive residues first.
    assert [group for group, _ in options] == ['Lys (K)', 'Glu (E)']
    assert options[0][1] == [('Lys-Arg', '→ Arg (R) · 2')]
    assert options[1][1] == [('Glu-Gly', '→ Gly (G) · 1'), ('Glu-Ter', '→ Stop · 1')]
    keys = [key for _, group in options for key, _ in group]
    assert sorted(keys) == sorted(table.matrix_counts)


def test_substitution_parsing_is_case_insensitive_and_bounded():
    assert coverage.parse_substitution('glu-GLY') == 'Glu-Gly'
    assert coverage.parse_substitution('Xyz-Gly') is None
    assert coverage.parse_substitution('Glu-Gly-Ala') is None
    assert coverage.parse_substitution('') is None


def _designs_table(rows, geometry=PLUS):
    """A ResultTable with its coverage attached the way the app does."""
    table = frame(rows)
    table.coverage = coverage.Coverage(table, geometry)
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
