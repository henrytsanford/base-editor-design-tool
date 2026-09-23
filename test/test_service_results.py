"""Parsing, filtering, sorting and paging a stored result.

Offline: the golden MAP2K1 designs file is a real engine output with 560 guides,
both strands, ClinVar matches and multi-edit guides, which is everything the table
has to cope with. No bundle and no HTTP.
"""
import threading
import time

import pytest

from bedesign import DESIGN_COLUMNS
from service.jobs import _tsv_gz
from service.params import TableView
from service.results import (ANY_MATCH, ROWS_PER_PAGE, ResultCache, ResultTable,
                             UnknownFilterValue)

GOLDEN = 'test/golden/map2k1/sgrna_designs.txt.gz'


@pytest.fixture(scope='module')
def raw():
    with open(GOLDEN, 'rb') as fh:
        return fh.read()


@pytest.fixture(scope='module')
def table(raw):
    return ResultTable(raw, DESIGN_COLUMNS)


def tsv_gz(rows):
    """A designs file with the real header and the given rows."""
    padded = [list(row) + [''] * (len(DESIGN_COLUMNS) - len(row)) for row in rows]
    return _tsv_gz(DESIGN_COLUMNS, padded)


def test_every_cell_survives_parsing_as_written(table):
    """dtype=str and na_filter=False, or the table would disagree with the download."""
    assert table.columns == DESIGN_COLUMNS
    assert table.total == 560
    page = table.select(TableView())
    assert all(isinstance(cell, str) for cell in page.rows[0])
    # An unedited guide has empty trailing fields, not 'nan' and not 0.
    assert page.rows[0][DESIGN_COLUMNS.index('Mutation category')] == ''


def test_a_multi_edit_guide_matches_each_of_its_categories(table):
    """'Missense;Silent;' is two edits, so it belongs to both filters."""
    missense = table.select(TableView(mutation='Missense')).matched
    silent = table.select(TableView(mutation='Silent')).matched
    both = table.select(TableView(mutation='Missense')).rows
    assert missense and silent
    index = DESIGN_COLUMNS.index('Mutation category')
    assert any('Silent' in row[index] for row in both), (
        'a guide categorised Missense;Silent; must appear under Missense')


def test_a_significance_filter_matches_whole_tokens_only(table):
    """'Likely benign' and 'Benign/Likely benign' are different classifications."""
    facets = dict(table.facets()['significance'])
    assert 'Likely benign' in facets
    likely = table.select(TableView(significance='Likely benign'))
    index = DESIGN_COLUMNS.index('Clinical significance')
    for row in likely.rows:
        tokens = row[index].split(';')
        assert 'Likely benign' in tokens
    assert likely.matched == facets['Likely benign']


def test_the_no_match_sentinel_is_never_offered_as_a_classification(table):
    """'None' means no ClinVar variant matched; it is not something to filter on."""
    assert 'None' not in dict(table.facets()['significance'])


def test_any_match_finds_the_rows_with_a_real_classification(table):
    facets = table.facets()
    any_hit = table.select(TableView(significance=ANY_MATCH))
    assert any_hit.matched == facets['any_match']
    index = DESIGN_COLUMNS.index('Clinical significance')
    for row in any_hit.rows:
        assert set(row[index].split(';')) - {'None', ''}


def test_a_filter_this_result_cannot_satisfy_is_refused(table):
    """Refused rather than silently empty, the way an unknown parameter is."""
    with pytest.raises(UnknownFilterValue):
        table.select(TableView(mutation='Nonexistent'))


def test_the_flag_filters_hide_only_the_flagged_rows(table):
    """The flag columns hold 'yes' or nothing -- never 'no'."""
    index = DESIGN_COLUMNS.index('4T flag')
    hidden = table.select(TableView(hide_4t=True))
    assert hidden.matched < table.total
    for row in hidden.rows:
        assert row[index] != 'yes'


def test_filters_combine(table):
    both = table.select(TableView(mutation='Missense', strand='sense')).matched
    one = table.select(TableView(mutation='Missense')).matched
    assert 0 < both <= one


def test_edit_counts_sort_as_numbers_not_text(table):
    """Every cell is a string, so without coercion '10' would sort before '9'."""
    index = DESIGN_COLUMNS.index('# edits')
    page = table.select(TableView(sort='# edits', dir='desc'))
    counts = [int(row[index]) for row in page.rows]
    assert counts == sorted(counts, reverse=True)
    assert counts[0] > 1


def test_a_blank_numeric_cell_sorts_last_rather_than_as_zero(table):
    rows = [['gA'] + [''] * 23, ['gB'] + [''] * 23]
    rows[0][DESIGN_COLUMNS.index('# edits')] = '2'
    small = ResultTable(tsv_gz(rows), DESIGN_COLUMNS)
    page = small.select(TableView(sort='# edits', dir='asc'))
    assert [row[0] for row in page.rows] == ['gA', 'gB']


def test_a_page_holds_fifty_rows_and_the_next_page_does_not_repeat_them(table):
    first = table.select(TableView(page=1))
    second = table.select(TableView(page=2))
    assert len(first.rows) == ROWS_PER_PAGE
    assert first.pages == 12 and second.page == 2
    assert not {row[0] for row in first.rows} & {row[0] for row in second.rows}


def test_a_page_past_the_end_lands_on_the_last_one(table):
    """A link into a result that is now filtered down should still show rows."""
    page = table.select(TableView(page=9999))
    assert page.page == page.pages and page.rows


def test_an_empty_result_still_renders(table):
    empty = ResultTable(tsv_gz([]), DESIGN_COLUMNS)
    page = empty.select(TableView())
    assert empty.columns == DESIGN_COLUMNS
    assert page.rows == [] and page.matched == 0 and page.pages == 1


def test_the_cache_parses_once_per_key(raw):
    cache = ResultCache()
    calls = []

    def load():
        calls.append(1)
        return ResultTable(raw, DESIGN_COLUMNS)

    first = cache.get('a' * 64, load)
    assert cache.get('a' * 64, load) is first
    assert len(calls) == 1


def test_simultaneous_misses_on_one_key_share_one_parse(raw):
    """Parsing separately would multiply a large result's memory by the number of
    requests waiting for it."""
    cache = ResultCache()
    calls = []
    release = threading.Event()

    def load():
        calls.append(1)
        release.wait(5)
        return ResultTable(raw, DESIGN_COLUMNS)

    got = []
    threads = [threading.Thread(target=lambda: got.append(cache.get('k', load)))
               for _ in range(8)]
    for thread in threads:
        thread.start()
    time.sleep(0.2)  # let every thread reach the cache before the parse finishes
    release.set()
    for thread in threads:
        thread.join(5)
    assert len(calls) == 1
    assert len(got) == 8 and all(table is got[0] for table in got)


def test_a_failed_parse_reaches_every_waiter_and_is_retried(raw):
    cache = ResultCache()
    release = threading.Event()

    def fails():
        release.wait(5)
        raise KeyError('evicted')

    errors = []

    def view():
        try:
            cache.get('k', fails)
        except KeyError as e:
            errors.append(e)

    threads = [threading.Thread(target=view) for _ in range(3)]
    for thread in threads:
        thread.start()
    time.sleep(0.2)
    release.set()
    for thread in threads:
        thread.join(5)
    assert len(errors) == 3
    # Nothing was cached, and nothing is left waiting: the next view parses afresh.
    assert cache.get('k', lambda: ResultTable(raw, DESIGN_COLUMNS)).total == 560


def test_the_cache_evicts_by_rows_not_by_entries(raw):
    """A frame costs ~1.3 KB a row, so entry count alone would not bound memory."""
    cache = ResultCache(max_rows=600, max_frames=8)
    for i in range(3):
        cache.get('key%d' % i, lambda: ResultTable(raw, DESIGN_COLUMNS))
    stats = cache.stats()
    assert stats['frames'] == 1 and stats['rows'] == 560


def test_the_cache_keeps_a_frame_larger_than_its_whole_budget(raw):
    """Otherwise the biggest genes, which need the cache most, never get it."""
    cache = ResultCache(max_rows=1, max_frames=8)
    cache.get('big', lambda: ResultTable(raw, DESIGN_COLUMNS))
    assert cache.stats()['frames'] == 1


def test_the_cache_drops_the_least_recently_used_frame(raw):
    cache = ResultCache(max_rows=10000, max_frames=2)
    for key in ('a', 'b'):
        cache.get(key, lambda: ResultTable(raw, DESIGN_COLUMNS))
    cache.get('a', lambda: ResultTable(raw, DESIGN_COLUMNS))
    reloaded = []
    cache.get('c', lambda: reloaded.append(1) or ResultTable(raw, DESIGN_COLUMNS))
    cache.get('a', lambda: reloaded.append(1) or ResultTable(raw, DESIGN_COLUMNS))
    assert len(reloaded) == 1, 'a was used most recently and should have survived'
