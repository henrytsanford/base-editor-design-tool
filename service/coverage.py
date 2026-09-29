"""Where the guides are, and what they can do to the protein.

Two pictures of one result, computed together because they read the same two columns
and answer to the same two filters.

The *map* lays the transcript out left to right with its introns and UTRs shortened
to the buffer the design actually tiled into, stacks the guides into a fixed number of
bins, and colours each bin by the worst thing its guides do. Binning is what keeps the
drawing the same size for TTN as for MAP2K1: the shape count follows the width of the
figure, never the guide count.

The *matrix* is the substitution table the base-editing papers draw -- initial residue
down, resulting residue across, grouped by chemical class. A filled square means this
editor can make that change somewhere in this transcript. Which squares fill is a
property of the genetic code and the editor rather than of the gene, so the same
editor fills the same squares in every transcript; only the counts and the positions
differ, which is what the map above it shows.

Both panels hand back plain data. The SVG itself is written in `_gene_map.html` and
`_matrix.html`, so the markup lives with the other markup and what is tested here is
arithmetic.
"""
import re
from collections import defaultdict

import numpy as np
import pandas as pd

from bedesign.engine import get_aa_map

from .results import (AA_COLUMN, CLASS_LABELS, CLASS_RANK, CLASSES, POSITION_COLUMN,
                      SIGNIFICANCE_COLUMN, STRAND_COLUMN, split_tokens)

# How many bins the transcript is divided into. The map's shape count is bounded by
# this and by the exon count, so a 35,991-residue gene draws about as much as a 393.
BINS = 300

# Drawing geometry, in the SVG's own user units. The template scales it to whatever
# width the page gives it; the map spans the page, so it is kept flat.
WIDTH = 1000
HEIGHT = 150
MARGIN = 10           # either side of the gene map
MIDLINE = 75
BAR_MAX = 46          # tallest a bin's stack may be drawn
TRACK = 13            # from the midline to where each strand's stacks start
CDS_HEIGHT = 16       # a coding exon block
UTR_HEIGHT = 8        # an untranslated one, drawn thinner the way a gene model is
LABEL_GAP = 22        # how far apart exon numbers must be to both be drawn
# What the unreachable middle of a long intron or UTR is squeezed into, in bases of
# drawing: enough to read as a break without costing the exons any width.
STUB = 24

# The consequence class the map does not stack: a guide with no edit in its window.
UNDRAWN = 'none'

# ClinVar classifications that make a guide worth flagging on the map.
PATHOGENIC = frozenset(['Pathogenic', 'Likely pathogenic',
                        'Pathogenic/Likely pathogenic'])

# A residue change as the engine writes it, e.g. 'Glu27Gly'. Splice edits are written
# 'Exon2:+1' instead and do not match, which is how they stay out of the matrix.
AA_EDIT = re.compile(r'^([A-Za-z]{3})(\d+)([A-Za-z]{3})$')

# The substitution matrix' axes, in the order the base-editing papers group them.
RESIDUE_GROUPS = (
    ('Pos.', ('Lys', 'Arg', 'His')),
    ('Neg.', ('Asp', 'Glu')),
    ('Polar', ('Ser', 'Thr', 'Asn', 'Gln')),
    ('Other', ('Cys', 'Gly', 'Pro')),
    ('Nonpolar', ('Ala', 'Val', 'Ile', 'Leu', 'Met')),
    ('Arom.', ('Phe', 'Tyr', 'Trp')),
    # Unnamed: its one residue is labelled 'Stop' already.
    ('', ('Ter',)),
)
RESIDUES = tuple(aa for _, group in RESIDUE_GROUPS for aa in group)
RESIDUE_INDEX = {aa: i for i, aa in enumerate(RESIDUES)}
ONE_LETTER = dict(get_aa_map(), Ter='Stop')

# Matrix geometry. Every label sits left of or below the grid, so the top margin is
# only enough to keep the frame's stroke inside the drawing.
CELL = 18
MATRIX_LEFT = 104
MATRIX_TOP = 6
# A square is coloured by what the change does, in the map's consequence hues, so one
# colour means one thing on both charts. How many guides make it is the fill's
# opacity, from MIN_OPACITY for the fewest to solid for the most: the palest square
# still has to read as filled beside an empty one.
MIN_OPACITY = 0.25
# How many steps the legend draws between the two ends of the opacity scale.
LEGEND_STEPS = 5


def square_class(source, target):
    """The consequence class of changing `source` to `target`, as the map names it."""
    if source == target:
        return 'sil'
    return 'lof' if target == 'Ter' else 'mis'


def plural(count, word):
    """'1 guide', '2 guides'."""
    return '%d %s%s' % (count, word, '' if count == 1 else 's')


def opacity(fraction):
    """The fill opacity at `fraction` (0 to 1) along the count scale."""
    fraction = min(1.0, max(0.0, fraction))
    return round(MIN_OPACITY + (1 - MIN_OPACITY) * fraction, 3)

# A substitution as a URL says it: 'Glu-Gly'. Both halves must name a residue, which
# is what stops the parameter reaching anything but a dictionary lookup.
SUBSTITUTION = re.compile(r'^([A-Za-z]{3})-([A-Za-z]{3})$')


def substitution_key(source, target):
    """The canonical name for one square, used in URLs and as a dict key."""
    return '%s-%s' % (source, target)


def parse_substitution(value):
    """`value` as its canonical key, e.g. 'glu-gly' as 'Glu-Gly', or None if it
    names no square."""
    match = SUBSTITUTION.match(value or '')
    if not match:
        return None
    source, target = match.group(1).title(), match.group(2).title()
    if source not in RESIDUE_INDEX or target not in RESIDUE_INDEX:
        return None
    return substitution_key(source, target)


class Geometry(object):
    """A transcript's exons laid out with its introns and UTRs shortened.

    Coordinates are *oriented*: genomic position times strand, so that a position
    increasing means moving 5' to 3' whichever strand the gene is on, and one set of
    arithmetic serves both. `buffer` is the design's own intron buffer, so the drawing
    shows exactly the non-coding sequence the run was allowed to tile into and no more.
    """

    def __init__(self, exons, cds, strand, buffer=30):
        self.strand = 1 if strand >= 0 else -1
        buffer = max(buffer, 1)
        self.exons = [tuple(sorted((a * self.strand, b * self.strand)))
                      for a, b in exons]
        self.exons.sort()
        self.cds = sorted(tuple(sorted((a * self.strand, b * self.strand)))
                          for a, b in (cds or []))
        self._lows = np.array([low for low, _ in self.exons], dtype=np.float64)
        self._highs = np.array([high for _, high in self.exons], dtype=np.float64)
        # Each exon's (low, high, is_coding) pieces, so UTR draws thinner.
        self.parts = self._split_parts()
        knots = self._knots(buffer)
        self._genomic = np.array([g for g, _ in knots], dtype=np.float64)
        self._drawn = np.array([x for _, x in knots], dtype=np.float64)
        self.width = max(knots[-1][1], 1)

    def _knots(self, buffer):
        """The drawing as (oriented position, drawing position) breakpoints, with
        straight lines between them.

        Coding sequence runs base for base. Anything else -- an intron, a UTR, the
        flank past either end -- keeps `buffer` bases at each edge, the most a guide
        can reach, and has the rest squeezed into a fixed stub, so a 5 kb UTR or
        intron draws no wider than a 100 bp one and cannot crowd out the exons.
        Each base spans [p, p + 1), so an exon ends at high + 1.
        """
        knots = [(self.exons[0][0] - buffer, 0.0)]

        def run(end):
            start, x = knots[-1]
            if end > start:
                knots.append((end, x + end - start))

        def squeeze(end):
            start, x = knots[-1]
            if end - start > 2 * buffer + STUB:
                knots.append((start + buffer, x + buffer))
                knots.append((end - buffer, x + buffer + STUB))
            run(end)

        for index, (low, _) in enumerate(self.exons):
            (squeeze if index else run)(low)
            for _, end, coding in self.parts[index]:
                (run if coding else squeeze)(end + 1)
        run(self.exons[-1][1] + 1 + buffer)
        return knots

    def x(self, positions):
        """Oriented coordinates as positions along the drawing, clamped to it."""
        return np.interp(positions, self._genomic, self._drawn)

    def place(self, positions):
        """Oriented coordinates as (positions along the drawing, nearest exon index).

        Vectorised because it is called once per result with every guide's position.
        """
        positions = np.asarray(positions, dtype=np.float64)
        count = len(self.exons)
        index = np.clip(np.searchsorted(self._highs, positions, side='left'),
                        0, count - 1)
        previous = np.clip(index - 1, 0, count - 1)
        before = self._lows[index] - positions
        after = positions - self._highs[previous]

        # The nearer of the two exons it sits between; a tie goes to the earlier one.
        # Past the last exon `before` is negative, which counts as inside it.
        distance = np.maximum(before, 0.0)
        exon = np.where(distance < np.abs(after), index, previous)
        return self.x(positions), exon.astype(np.int32)

    def _split_parts(self):
        """Every exon's pieces, in one merge pass: exons and CDS spans are both
        sorted, so each span is visited only by the exons it overlaps."""
        split, first = [], 0
        for low, high in self.exons:
            while first < len(self.cds) and self.cds[first][1] < low:
                first += 1
            pieces, cursor, span = [], low, first
            while span < len(self.cds) and self.cds[span][0] <= high:
                start = max(low, self.cds[span][0])
                end = min(high, self.cds[span][1])
                if start > cursor:
                    pieces.append((cursor, start - 1, False))
                pieces.append((start, end, True))
                cursor = end + 1
                span += 1
            if cursor <= high:
                pieces.append((cursor, high, False))
            split.append(pieces)
        return split


class Coverage(object):
    """The two panels, and the masks that let the table answer to them."""

    def __init__(self, table, geometry):
        """Built over a ResultTable, whose own consequence classes are reused so
        the map and the filter cannot count a guide differently."""
        self.geometry = geometry
        self.total = table.total
        self._scan(table.frame, table.row_classes)
        self.matrix_counts = {key: len(rows)
                              for key, rows in self._substitution_rows.items()}
        n_exons = len(geometry.exons)
        drawn = self._row_exons >= 0
        per_bin = np.bincount(self._row_bins[drawn] * n_exons + self._row_exons[drawn],
                              minlength=BINS * n_exons).reshape(BINS, n_exons)
        self.exon_guides = per_bin.sum(axis=0)
        # The exon each bin's click filters to, numbered from 1 as the exon filter
        # takes it: the one most of its guides sit in, since a bin at a splice
        # junction can straddle two. 0 for an empty bin.
        self.bin_exons = (np.where(per_bin.any(axis=1), per_bin.argmax(axis=1), -1)
                          + 1).tolist()
        # How many of each substitution's guides every exon holds, so either chart
        # selection can mark the pairs with the other that would match nothing.
        self._substitution_exons = {}
        for key, rows in self._substitution_rows.items():
            exons_of = self._row_exons[rows]
            self._substitution_exons[key] = np.bincount(exons_of[exons_of >= 0],
                                                        minlength=n_exons)
        self.peak = max(int(self._bins.sum(axis=2).max()), 1)

    def _scan(self, frame, worst):
        """One pass over the frame, filling both panels.

        Read column by column as numpy arrays rather than row by row: on TTN's 25,662
        guides the difference is the whole cost of the feature.
        """
        positions = pd.to_numeric(frame[POSITION_COLUMN], errors='coerce').to_numpy()
        senses = (frame[STRAND_COLUMN].to_numpy() == 'sense')
        amino = frame[AA_COLUMN].to_numpy()
        significance = frame[SIGNIFICANCE_COLUMN].to_numpy()

        # Placed in one pass; a row whose position did not parse is dropped from the
        # drawing but still counted in the table, so the two never disagree about how
        # many guides there are. Per row, the exon (-1 when undrawn) and the bin are
        # kept so a selected substitution can be marked without walking the frame.
        drawable = np.isfinite(positions)
        oriented = np.where(drawable, positions, 0.0) * self.geometry.strand
        places, exons = self.geometry.place(oriented)
        self._row_bins = np.clip((places / self.geometry.width * BINS).astype(np.int32),
                                 0, BINS - 1)
        self._row_exons = np.where(drawable, exons, -1).astype(np.int32)

        # Kept per row so a selected substitution's guides can be lifted out of the
        # stacks they were counted into.
        self._row_sides = (~senses).astype(np.int32)
        self._row_classes = worst
        # Counts per (bin, strand side, class), sense above the line as side 0. A
        # guide that edits nothing has no consequence to draw, so it is left off the
        # stacks; the legend still counts it and the table still lists it.
        stacked = drawable & (worst != CLASS_RANK[UNDRAWN])
        cell = ((self._row_bins * 2 + self._row_sides) * len(CLASSES)
                + worst)[stacked]
        self._bins = np.bincount(cell, minlength=BINS * 2 * len(CLASSES)).reshape(
            BINS, 2, len(CLASSES))

        substitution_rows, sites = defaultdict(list), defaultdict(set)
        self.matrix_pathogenic = defaultdict(int)
        for row in np.flatnonzero(drawable):
            edits = split_tokens(amino[row])
            marks = split_tokens(significance[row])
            for position_in_list, edit in enumerate(edits):
                match = AA_EDIT.match(edit)
                if not match:
                    continue
                source, site, target = match.group(1), int(match.group(2)), match.group(3)
                if source not in RESIDUE_INDEX or target not in RESIDUE_INDEX:
                    continue
                key = substitution_key(source, target)
                sites[key].add(site)
                substitution_rows[key].append(row)
                mark = marks[position_in_list] if position_in_list < len(marks) else ''
                if mark in PATHOGENIC:
                    self.matrix_pathogenic[key] += 1
        # A guide that makes the same substitution at two residues is one guide, and
        # the square has to report the number the table will show when it is clicked.
        self._substitution_rows = {key: np.unique(np.asarray(rows, dtype=np.int32))
                                   for key, rows in substitution_rows.items()}
        self.matrix_sites = {key: len(found) for key, found in sites.items()}

    # ---- filters -------------------------------------------------------------

    def exon_mask(self, number):
        """Rows whose guide sits in exon `number`, counting from 1."""
        return self._row_exons == number - 1

    def substitution_mask(self, key):
        """Rows making the substitution `key`, e.g. 'Glu-Gly'."""
        mask = np.zeros(self.total, dtype=bool)
        mask[self._substitution_rows[key]] = True
        return mask

    def has_exon(self, number):
        return 1 <= number <= len(self.geometry.exons)

    def has_substitution(self, key):
        return key in self.matrix_counts

    # ---- what the template draws ---------------------------------------------

    def map_view(self, selected_exon=0, selected_substitution='', consequence=''):
        """The gene map: exon blocks and binned guide stacks.

        With a `consequence` selected, every segment of another class is marked `dim`,
        and a bin holding none of that class is marked `dim` whole.

        With a substitution selected, the guides making it are drawn as their own
        segment next to the track, in the colour of the square that was clicked, and
        every other segment is marked `dim`. They take the square's colour rather
        than their own because a guide is otherwise drawn as its worst edit: one
        making a silent change beside a missense one would light up blue when the
        green square was clicked. A bin holding none of them is marked `dim` whole.
        """
        geometry = self.geometry
        scale = (WIDTH - 2 * MARGIN) / geometry.width
        bin_width = (WIDTH - 2 * MARGIN) / BINS
        # With a substitution selected, how many of its guides each exon holds: an
        # exon holding none would filter the table to nothing, so it is marked
        # `empty` and the page offers no link off it.
        exon_picked = self._substitution_exons.get(selected_substitution)
        exons = []
        last_label = -LABEL_GAP
        for index, (low, high) in enumerate(geometry.exons):
            pieces = geometry.parts[index]
            blocks = []
            for start, end, coding in pieces:
                x0 = MARGIN + geometry.x(start) * scale
                x1 = MARGIN + geometry.x(end + 1) * scale
                height = CDS_HEIGHT if coding else UTR_HEIGHT
                blocks.append({'x': round(x0, 2),
                               'width': round(max(0.8, x1 - x0), 2),
                               'y': round(MIDLINE - height / 2.0, 2),
                               'height': height})
            # The number sits on the coding part when there is one: white on the
            # thin UTR block it would be cut off or unreadable.
            coding_pieces = [piece for piece in pieces if piece[2]]
            if coding_pieces:
                first, last = coding_pieces[0][0], coding_pieces[-1][1]
            else:
                first, last = low, high
            centre = MARGIN + (geometry.x(first) + geometry.x(last + 1)) / 2.0 * scale
            label = ''
            if centre - last_label >= LABEL_GAP:
                label, last_label = str(index + 1), centre
            count = int(self.exon_guides[index])
            empty = exon_picked is not None and not exon_picked[index]
            exons.append({
                'number': index + 1,
                'blocks': blocks,
                'label': label,
                'label_x': round(centre, 2),
                'selected': selected_exon == index + 1,
                'empty': empty,
                'title': 'Exon %d · %d bp%s · %s%s' % (
                    index + 1, high - low + 1,
                    '' if coding_pieces else ' (untranslated)',
                    plural(count, 'guide'),
                    ' · none make the selected change' if empty else ''),
            })

        # The selected guides per (bin, side), and the same guides per (bin, side,
        # class) so they can be taken out of the counts they were stacked under.
        picked = removed = None
        if selected_substitution in self._substitution_rows:
            rows = self._substitution_rows[selected_substitution]
            # With a consequence selected too, only the guides both filters keep.
            if consequence:
                rows = rows[self._row_classes[rows] == CLASS_RANK[consequence]]
            picked_class = square_class(*selected_substitution.split('-'))
            picked_letters = self.selection(selected_substitution)['letters']
            where = (self._row_bins[rows], self._row_sides[rows])
            picked = np.zeros((BINS, 2), dtype=np.int64)
            np.add.at(picked, where, 1)
            removed = np.zeros_like(self._bins)
            np.add.at(removed, where + (self._row_classes[rows],), 1)
        bar_width = round(max(0.6, bin_width - 0.6), 2)
        bins = []
        for index, sides in enumerate(self._bins.tolist()):
            both = [up + down for up, down in zip(*sides)]
            total = sum(both)
            if not total:
                continue
            x = MARGIN + index * bin_width
            bar_x = round(x + 0.3, 2)
            segments = []
            for side, counts in enumerate(sides):
                if picked is None:
                    stack = [(name, value, bool(consequence) and name != consequence)
                             for name, value in zip(CLASSES, counts)]
                else:
                    stack = [(picked_class, int(picked[index, side]), False)] + [
                        (name, value - int(removed[index, side, rank]), True)
                        for rank, (name, value) in enumerate(zip(CLASSES, counts))]
                edge = MIDLINE - TRACK if side == 0 else MIDLINE + TRACK
                for name, value, faded in stack:
                    if not value:
                        continue
                    height = value / self.peak * BAR_MAX
                    y = edge - height if side == 0 else edge
                    segments.append({'y': round(y, 2),
                                     'height': round(height, 2),
                                     'cls': name,
                                     'dim': faded})
                    edge = edge - height if side == 0 else edge + height
            hits = int(picked[index].sum()) if picked is not None else 0
            bins.append({
                'x': round(x, 2),
                'bar_x': bar_x,
                'segments': segments,
                'dim': all(segment['dim'] for segment in segments),
                'exon': self.bin_exons[index],
                'title': '%s: %s%s' % (
                    plural(total, 'guide'),
                    ', '.join('%d %s' % (n, CLASS_LABELS[c])
                              for c, n in zip(CLASSES, both) if n),
                    ' · %d make %s' % (hits, picked_letters) if hits else ''),
            })

        return {'exons': exons, 'bins': bins,
                'width': WIDTH, 'height': HEIGHT,
                'bin_width': round(bin_width, 2), 'bar_width': bar_width,
                # The band the stacks can occupy: the hover targets span exactly
                # this, whatever the drawing's height.
                'band_top': MIDLINE - TRACK - BAR_MAX,
                'band_height': 2 * (TRACK + BAR_MAX),
                'midline': MIDLINE, 'margin': MARGIN}

    def matrix_view(self, selected='', exon=0):
        """The substitution matrix: 441 squares, whatever the gene.

        With an exon selected, a square none of whose guides sit in that exon is
        marked `empty`: the two selections combine, and that pair matches nothing.
        """
        # The scale spans the squares that change the residue. Silent ones sit on it
        # too but are left out of its range, so a large synonymous count does not
        # wash out the rest; past the top they are simply drawn solid.
        shaded = [count for key, count in self.matrix_counts.items()
                  if square_class(*key.split('-')) != 'sil']
        low, high = min(shaded + [1]), max(shaded + [1])
        cells, present = [], set()
        for down, source in enumerate(RESIDUES):
            for across, target in enumerate(RESIDUES):
                key = substitution_key(source, target)
                count = self.matrix_counts.get(key, 0)
                shade = None
                empty = bool(count and exon) and not self._substitution_exons[key][exon - 1]
                if count:
                    cls = square_class(source, target)
                    present.add(cls)
                    shade = opacity((count - low) / float(high - low)
                                    if high > low else 1.0)
                    sites = self.matrix_sites[key]
                    pathogenic = self.matrix_pathogenic.get(key, 0)
                    title = '%s to %s · %s at %s%s%s' % (
                        ONE_LETTER[source], ONE_LETTER[target],
                        plural(count, 'guide'), plural(sites, 'site'),
                        ' · %d recreate a pathogenic variant' % pathogenic
                        if pathogenic else '',
                        ' · none in exon %d' % exon if empty else '')
                else:
                    cls, title = 'off', '%s to %s · not reachable with this editor' % (
                        ONE_LETTER[source], ONE_LETTER[target])
                cells.append({
                    'key': key,
                    'x': MATRIX_LEFT + across * CELL,
                    'y': MATRIX_TOP + down * CELL,
                    'count': count,
                    'cls': cls,
                    'opacity': shade,
                    'selected': selected == key,
                    'empty': empty,
                    'title': title,
                })

        rules, groups = [], []
        edge = 0
        for name, members in RESIDUE_GROUPS:
            start = edge
            edge += len(members)
            rules.append(edge)
            if name:
                groups.append({'name': name,
                               'centre': round(start + len(members) / 2.0, 3)})
        labels = [ONE_LETTER[residue] for residue in RESIDUES]

        size = len(RESIDUES) * CELL
        return {
            'cells': cells, 'rules': rules, 'labels': labels, 'groups': groups,
            'left': MATRIX_LEFT, 'top': MATRIX_TOP, 'cell': CELL, 'size': size,
            'width': MATRIX_LEFT + size + 14, 'height': MATRIX_TOP + size + 72,
            # The legend: the classes this result fills, in the map's order and
            # words, and the opacity scale as a row of steps.
            'classes': [{'cls': name, 'label': CLASS_LABELS[name]}
                        for name in CLASSES if name in present],
            'steps': [opacity(i / float(LEGEND_STEPS - 1))
                      for i in range(LEGEND_STEPS)],
            'low': low, 'high': high, 'shaded': bool(shaded),
        }

    def selection(self, key):
        """A one-line summary of the chosen square, for the status line."""
        source, target = key.split('-')
        return {
            'letters': '%s → %s' % (ONE_LETTER[source], ONE_LETTER[target]),
            # The filter's chip, in the three-letter names the table's edits use.
            'label': '%s → %s' % (source, 'Stop' if target == 'Ter' else target),
            'sites': self.matrix_sites[key],
            'pathogenic': self.matrix_pathogenic.get(key, 0),
            'silent': source == target,
        }
