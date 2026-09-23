"""Parsed design results, and the bounded cache of them.

Design doc 5.3 measured the weakness this fixes: every table view re-gunzipped and
re-parsed the whole TSV -- 17 ms for TTN at 12,842 guides, ~35 ms at `edit=all`'s
25,662 -- and with no cache, concurrent views serialise on the GIL. Parsing once and
keeping the frame turns every later view into a mask lookup and a slice.

Two columns are not scalars. `Mutation category` and `Clinical significance` each hold
one ';'-joined entry per edit in the window, so a guide with two edits reads
'Missense;Silent;'. Filtering therefore matches a *token*, never the whole cell, and
the vocabularies are read off the data rather than written down here: the mutation
categories are a closed set the engine owns, but clinical significance is whatever
ClinVar's classification column says, which changes without asking us.

Nothing here knows about HTTP. `select` takes the view object `params.parse_view_query`
built and gives back one page of rows.
"""
import io
import math
import threading
from collections import OrderedDict
from concurrent.futures import Future
from dataclasses import dataclass

import numpy as np
import pandas as pd

# Design doc 5: 50 rows per page.
ROWS_PER_PAGE = 50

MUTATION_COLUMN = 'Mutation category'
SIGNIFICANCE_COLUMN = 'Clinical significance'
# The two ';'-joined columns, the only ones filtered by token.
TOKEN_COLUMNS = (MUTATION_COLUMN, SIGNIFICANCE_COLUMN)

# What the engine writes into 'Clinical significance' for an edit that matched no
# ClinVar SNP. A sentinel rather than a classification, so it is never offered as a
# facet -- 'no match' is the absence of one.
NO_MATCH = 'None'
# Asks for rows carrying any real classification at all. ClinVar's vocabulary is open,
# but it is a controlled list of clinical terms, so this cannot collide with one.
ANY_MATCH = 'any-match'

# Filtered by an exact value rather than by a token. Two distinct values each, so
# masking them once per parse costs almost nothing and takes the elementwise string
# comparison off every filtered request.
VALUE_COLUMNS = ('Edit', 'sgRNA Strand', 'BsmBI flag', '4T flag')

# Sorted as numbers. Every cell is a string, so without this '10' sorts before '9'.
NUMERIC_COLUMNS = frozenset(['# edits', '#silent edits', 'sgrna genomic position'])

# The flag columns hold 'yes' or the empty string -- never 'no'.
FLAG_YES = 'yes'


class UnknownFilterValue(ValueError):
    """A filter naming a value this result does not contain.

    Open-vocabulary filters cannot be checked when the query string is parsed, only
    against the result they are applied to, so this is raised here and rendered as a
    400 rather than quietly returning no rows.
    """


@dataclass(frozen=True)
class Page:
    """One page of a filtered, sorted result."""
    rows: list
    matched: int
    page: int
    pages: int


def _masks_from(positions, length):
    masks = {}
    for value, rows in positions.items():
        mask = np.zeros(length, dtype=bool)
        mask[rows] = True
        masks[value] = mask
    return masks


def _token_masks(values):
    """Boolean masks, one per token appearing in a ';'-joined column.

    Built from the split cells rather than a substring test: 'Benign',
    'Likely benign' and 'Benign/Likely benign' are three distinct ClinVar
    classifications, and a substring test would make the first match all three.
    """
    positions = {}
    for row, cell in enumerate(values):
        if not cell:
            continue
        for token in cell.split(';'):
            if token:
                positions.setdefault(token, []).append(row)
    return _masks_from(positions, len(values))


def _value_masks(values):
    """Boolean masks, one per distinct value in a scalar column."""
    positions = {}
    for row, cell in enumerate(values):
        positions.setdefault(cell, []).append(row)
    return _masks_from(positions, len(values))


class ResultTable(object):
    """One parsed designs.tsv.gz, with its filter masks precomputed."""

    def __init__(self, raw, columns):
        try:
            self.frame = pd.read_csv(
                io.BytesIO(raw), sep='\t', compression='gzip',
                # Every cell stays the exact text the engine wrote, so the table shows
                # what the download contains. Without this pandas turns '' into NaN
                # and '00123' into a number.
                dtype=str, keep_default_na=False, na_filter=False)
        except pd.errors.EmptyDataError:
            self.frame = pd.DataFrame(columns=list(columns))
        self.columns = list(self.frame.columns)
        self.total = len(self.frame)
        self._masks = {
            column: _token_masks(self.frame[column].to_numpy())
            for column in TOKEN_COLUMNS if column in self.frame.columns}
        self._masks.update(
            (column, _value_masks(self.frame[column].to_numpy()))
            for column in VALUE_COLUMNS if column in self.frame.columns)
        self._any_match = self._compute_any_match()

    def _compute_any_match(self):
        """Rows carrying at least one real ClinVar classification."""
        masks = self._masks.get(SIGNIFICANCE_COLUMN, {})
        found = np.zeros(self.total, dtype=bool)
        for token, mask in masks.items():
            if token != NO_MATCH:
                found |= mask
        return found

    def facets(self):
        """The filter vocabularies present in *this* result, with row counts.

        Read off the data, so a ClinVar classification we have never seen shows up in
        the dropdown instead of being silently dropped by a hardcoded list.
        """
        def counted(column, skip=()):
            masks = self._masks.get(column, {})
            return sorted((token, int(mask.sum()))
                          for token, mask in masks.items() if token not in skip)

        return {
            'mutation': counted(MUTATION_COLUMN),
            'significance': counted(SIGNIFICANCE_COLUMN, skip=(NO_MATCH,)),
            'any_match': int(self._any_match.sum()),
        }

    def _token_mask(self, column, token):
        masks = self._masks.get(column, {})
        if token not in masks:
            raise UnknownFilterValue(token)
        return masks[token]

    def _value_mask(self, column, value):
        """Rows whose cell is exactly `value`.

        A value this column does not contain is an empty mask, not an error: the
        vocabularies here are closed and already checked when the query was parsed,
        so 'no guide is on the antisense strand' is an answer rather than a mistake.
        """
        mask = self._masks.get(column, {}).get(value)
        return np.zeros(self.total, dtype=bool) if mask is None else mask

    def _filter_mask(self, view):
        mask = None

        def keep(other):
            return other if mask is None else (mask & other)

        if view.mutation:
            mask = keep(self._token_mask(MUTATION_COLUMN, view.mutation))
        if view.significance:
            if view.significance == ANY_MATCH:
                mask = keep(self._any_match)
            else:
                mask = keep(self._token_mask(SIGNIFICANCE_COLUMN, view.significance))
        if view.deaminase:
            mask = keep(self._value_mask('Edit', view.deaminase))
        if view.strand:
            mask = keep(self._value_mask('sgRNA Strand', view.strand))
        if view.hide_bsmbi:
            mask = keep(~self._value_mask('BsmBI flag', FLAG_YES))
        if view.hide_4t:
            mask = keep(~self._value_mask('4T flag', FLAG_YES))
        return mask

    def _ordered(self, positions, view):
        """`positions`, reordered by the sort this view asks for.

        Sorts the one column being sorted on rather than the frame: a page is 50 rows,
        so carrying every matched row across all 24 columns through a filter and a
        reindex is work that is thrown away.
        """
        if not view.sort or view.sort not in self.frame.columns or not positions.size:
            return positions
        keys = self.frame[view.sort].take(positions).reset_index(drop=True)
        if view.sort in NUMERIC_COLUMNS:
            # A blank cell is not zero -- it is 'this guide has no edits'. coerce puts
            # those at the end either way, rather than sorting them among the numbers.
            keys = pd.to_numeric(keys, errors='coerce')
        order = keys.sort_values(kind='mergesort', ascending=view.dir != 'desc',
                                 na_position='last')
        return positions[order.index.to_numpy()]

    def select(self, view):
        """One page of rows, after filtering and sorting.

        The page is clamped to the last one that exists, so a stale link deep into a
        result that is now filtered down still lands on rows rather than on nothing.
        """
        mask = self._filter_mask(view)
        positions = np.arange(self.total) if mask is None else np.flatnonzero(mask)
        matched = int(positions.size)
        pages = max(1, math.ceil(matched / ROWS_PER_PAGE))
        page = min(max(view.page, 1), pages)
        positions = self._ordered(positions, view)
        start = (page - 1) * ROWS_PER_PAGE
        window = self.frame.take(positions[start:start + ROWS_PER_PAGE])
        return Page(rows=window.values.tolist(), matched=matched, page=page,
                    pages=pages)


class ResultCache(object):
    """A bounded LRU of parsed results, shared by the ASGI thread pool.

    Bounded by rows rather than by entries: a frame costs about 1.3 KB per row, so
    TTN's 25,662 guides are ~34 MB and a count of entries would not be a memory bound
    at all. One frame is always admitted even if it alone exceeds the budget, since
    the alternative is that the biggest genes -- the ones that most need the cache --
    are the only ones that never get it. That frame is bounded all the same: the app
    never parses a result with more than `TABLE_MAX_ROWS` guides, and offers only
    its downloads instead.
    """

    def __init__(self, max_rows=100000, max_frames=8):
        self.max_rows = max_rows
        self.max_frames = max_frames
        self._lock = threading.Lock()
        self._tables = OrderedDict()
        self._rows = 0
        # key -> Future for a parse in progress, so misses on one key share it.
        self._loading = {}

    def get(self, key, load):
        """The parsed table for a key, calling `load()` only on a miss.

        `load` runs outside the lock, so parsing one result never makes views of
        another wait. Simultaneous misses on one key share a single `load()`: the
        first caller parses and the rest wait for its result. Parsing separately
        would multiply memory by the number of waiting requests, which for a large
        result is the difference between a slow page and running out of memory.

        If `load()` raises, every caller waiting on it gets the same exception and
        the next `get` tries again.
        """
        with self._lock:
            table = self._tables.get(key)
            if table is not None:
                self._tables.move_to_end(key)
                return table
            pending = self._loading.get(key)
            loading = pending is None
            if loading:
                pending = self._loading[key] = Future()
        if not loading:
            return pending.result()

        try:
            table = load()
        except BaseException as e:
            with self._lock:
                del self._loading[key]
            pending.set_exception(e)
            raise

        with self._lock:
            del self._loading[key]
            self._tables[key] = table
            self._rows += table.total
            self._tables.move_to_end(key)
            self._evict()
        pending.set_result(table)
        return table

    def _evict(self):
        """Called with the lock held. Never evicts the entry just inserted."""
        while len(self._tables) > 1 and (len(self._tables) > self.max_frames
                                         or self._rows > self.max_rows):
            _, evicted = self._tables.popitem(last=False)
            self._rows -= evicted.total

    def stats(self):
        with self._lock:
            return {'frames': len(self._tables), 'rows': self._rows}
