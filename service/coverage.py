"""Where the guides are, and what they can do to the protein.

Two pictures of one result, computed together because they read the same two columns
and answer to the same two filters.

The *map* lays the transcript out left to right with its introns shortened to the
buffer the design actually tiled into, stacks the guides into a fixed number of bins,
and colours each bin by the worst thing its guides do. Binning is what keeps the
drawing the same size for TTN as for MAP2K1: the shape count follows the width of the
figure, never the guide count.

The *matrix* is the substitution table the base-editing papers draw -- initial residue
down, resulting residue across, grouped by chemical class. A filled square means this
editor can make that change somewhere in this transcript. Which squares fill is a
property of the genetic code and the editor rather than of the gene, so the same
editor fills the same squares in every transcript; only the counts and the positions
differ, which is what the map beside it shows.

Both panels hand back plain data. The SVG itself is written in `_coverage.html`, so
the markup lives with the other markup and what is tested here is arithmetic.
"""
import re
from collections import defaultdict

import numpy as np
import pandas as pd

from bedesign.engine import get_aa_map

from .results import MUTATION_COLUMN, SIGNIFICANCE_COLUMN, split_tokens, token_masks

POSITION_COLUMN = 'sgrna genomic position'
STRAND_COLUMN = 'sgRNA Strand'
AA_COLUMN = 'Amino acid edits'

# How many bins the transcript is divided into. The map's shape count is bounded by
# this and by the exon count, so a 35,991-residue gene draws about as much as a 393.
BINS = 300

# Drawing geometry, in the SVG's own user units. The template scales it to whatever
# width the page gives it.
WIDTH = 1000
HEIGHT = 230
MARGIN = 10           # either side of the gene map
MIDLINE = 118
BAR_MAX = 78          # tallest a bin's stack may be drawn
CDS_HEIGHT = 16       # a coding exon block
UTR_HEIGHT = 8        # an untranslated one, drawn thinner the way a gene model is
LABEL_GAP = 22        # how far apart exon numbers must be to both be drawn

# What a guide does, worst first. A guide is drawn as the most severe edit it makes,
# so a guide that knocks out a splice site and also makes a silent change reads as
# the knockout.
CLASSES = ('lof', 'mis', 'sil', 'nc', 'none')
CLASS_RANK = {name: i for i, name in enumerate(CLASSES)}
CLASS_LABELS = {'lof': 'loss of function', 'mis': 'missense', 'sil': 'silent',
                'nc': 'non-coding edit', 'none': 'no edit in window'}
# The engine's categories, mapped onto them. Nonsense and a broken splice site are
# both a dead protein, so they share a colour.
CATEGORY_CLASS = {
    'Nonsense': 'lof', 'Splice-donor': 'lof', 'Splice-acceptor': 'lof',
    'Missense': 'mis', 'Silent': 'sil',
    'UTR': 'nc', 'Intron': 'nc', 'Flanking': 'nc',
}

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

# Matrix geometry. The square is inherently taller than the wide, flat gene map
# beside it, so the cell is kept small enough that the two panels balance.
CELL = 18
MATRIX_LEFT = 104
MATRIX_TOP = 52
# Count on the matrix, one hue light to dark. A square's fill is interpolated along
# these stops, and the legend draws the same stops as a gradient, so the two agree.
# Deliberately not one of the map's consequence hues: it encodes how many, not what.
RAMP = ('#cde2fb', '#86b6ef', '#3987e5', '#184f95')


def ramp_colour(fraction):
    """The colour at `fraction` (0 to 1) along RAMP."""
    fraction = min(1.0, max(0.0, fraction))
    span = fraction * (len(RAMP) - 1)
    index = min(int(span), len(RAMP) - 2)
    local = span - index
    low, high = RAMP[index], RAMP[index + 1]
    channels = (round(int(low[i:i + 2], 16) * (1 - local) + int(high[i:i + 2], 16) * local)
                for i in (1, 3, 5))
    return '#%02x%02x%02x' % tuple(channels)

# A substitution as a URL says it: 'Glu-Gly'. Both halves must name a residue, which
# is what stops the parameter reaching anything but a dictionary lookup.
SUBSTITUTION = re.compile(r'^([A-Za-z]{3})-([A-Za-z]{3})$')


def substitution_key(source, target):
    """The canonical name for one square, used in URLs and as a dict key."""
    return '%s-%s' % (source, target)


def parse_substitution(value):
    """`value` as a (source, target) pair, or None if it names no square."""
    match = SUBSTITUTION.match(value or '')
    if not match:
        return None
    source, target = match.group(1).title(), match.group(2).title()
    if source not in RESIDUE_INDEX or target not in RESIDUE_INDEX:
        return None
    return source, target


class Geometry(object):
    """A transcript's exons laid out with its introns shortened.

    Coordinates are *oriented*: genomic position times strand, so that a position
    increasing means moving 5' to 3' whichever strand the gene is on, and one set of
    arithmetic serves both. `buffer` is the design's own intron buffer, so the drawing
    shows exactly the intron the run was allowed to tile into and no more.
    """

    def __init__(self, exons, cds, strand, buffer=30):
        self.strand = 1 if strand >= 0 else -1
        buffer = max(buffer, 1)
        self.exons = [tuple(sorted((a * self.strand, b * self.strand)))
                      for a, b in exons]
        self.exons.sort()
        self.cds = sorted(tuple(sorted((a * self.strand, b * self.strand)))
                          for a, b in (cds or []))
        # An intron is drawn as the buffer either side plus a fixed stub, so a long
        # intron does not squeeze the exons into invisibility.
        self.gap = 2 * buffer + 24
        offsets = []
        x = buffer
        for low, high in self.exons:
            offsets.append(x)
            x += high - low + 1 + self.gap
        self.width = max(x - self.gap + buffer, 1)
        self.offsets = np.array(offsets, dtype=np.float64)
        self._lows = np.array([low for low, _ in self.exons], dtype=np.float64)
        self._highs = np.array([high for _, high in self.exons], dtype=np.float64)
        self._parts = self._split_parts()

    def place(self, positions):
        """Oriented coordinates as (positions along the drawing, nearest exon index).

        Vectorised because it is called once per result with every guide's position.
        """
        positions = np.asarray(positions, dtype=np.float64)
        count = len(self.exons)
        index = np.clip(np.searchsorted(self._highs, positions, side='left'),
                        0, count - 1)
        previous = np.clip(index - 1, 0, count - 1)
        low, high = self._lows[index], self._highs[index]
        previous_low, previous_high = self._lows[previous], self._highs[previous]
        half = self.gap / 2.0

        # Inside an exon the drawing runs base for base; past the last one it keeps
        # running, which is the same expression, so only the clamp is separate.
        inside = positions >= low
        within = self.offsets[index] + (positions - low)
        after = positions - previous_high
        before = low - positions
        upstream = (self.offsets[previous] + (previous_high - previous_low)
                    + np.minimum(after, half))
        downstream = self.offsets[index] - np.minimum(before, half)
        intron = np.where(after <= before, upstream, downstream)
        # Before the first exon there is no previous one to measure from.
        intron = np.where(index == 0,
                          np.maximum(0.0, self.offsets[0] - before), intron)
        layout = np.minimum(np.where(inside, within, intron), self.width)

        # The nearer of the two exons it sits between; a tie goes to the earlier one.
        distance = np.where(inside, 0.0, np.abs(before))
        exon = np.where(distance < np.abs(after), index, previous)
        return layout, exon.astype(np.int32)

    def parts(self, index):
        """(low, high, is_coding) pieces of one exon, so UTR draws thinner."""
        return self._parts[index]

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

    def __init__(self, frame, geometry, category_masks=None):
        """`category_masks` are the table's own masks for MUTATION_COLUMN, so the
        column is split once per result; built from `frame` when not given."""
        self.geometry = geometry
        self.total = len(frame)
        if category_masks is None:
            category_masks = token_masks(frame[MUTATION_COLUMN].to_numpy())
        self._scan(frame, category_masks)
        self.matrix_counts = {key: len(rows)
                              for key, rows in self._substitution_rows.items()}
        drawn = self._row_exons[self._row_exons >= 0]
        self.exon_guides = np.bincount(drawn, minlength=len(geometry.exons))
        self.peak = max(int(self._bins.sum(axis=2).max()), 1)

    def _scan(self, frame, category_masks):
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

        # Each row's worst class, as a rank into CLASSES: every category lowers the
        # rank of the rows carrying it to its own class's, if that is worse.
        worst = np.full(self.total, CLASS_RANK['none'], dtype=np.int32)
        for category, mask in category_masks.items():
            rank = CLASS_RANK[CATEGORY_CLASS.get(category, 'nc')]
            worst[mask] = np.minimum(worst[mask], rank)
        # Counts per (bin, strand side, class), sense above the line as side 0.
        cell = ((self._row_bins * 2 + (~senses).astype(np.int32)) * len(CLASSES)
                + worst)[drawable]
        self._bins = np.bincount(cell, minlength=BINS * 2 * len(CLASSES)).reshape(
            BINS, 2, len(CLASSES))
        self.class_counts = dict(zip(CLASSES, self._bins.sum(axis=(0, 1)).tolist()))

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

    def map_view(self, selected_exon=0, selected_substitution=''):
        """The gene map: exon blocks, binned guide stacks, and any highlight."""
        geometry = self.geometry
        scale = (WIDTH - 2 * MARGIN) / geometry.width
        bin_width = (WIDTH - 2 * MARGIN) / BINS
        exons = []
        last_label = -LABEL_GAP
        for index, (low, high) in enumerate(geometry.exons):
            pieces = geometry.parts(index)
            offset = geometry.offsets[index] - low
            blocks = []
            for start, end, coding in pieces:
                x0 = MARGIN + (offset + start) * scale
                x1 = MARGIN + (offset + end) * scale + scale
                height = CDS_HEIGHT if coding else UTR_HEIGHT
                blocks.append({'x': round(x0, 2),
                               'width': round(max(0.8, x1 - x0), 2),
                               'y': round(MIDLINE - height / 2.0, 2),
                               'height': height})
            centre = MARGIN + (geometry.offsets[index] + (high - low) / 2.0) * scale
            label = ''
            if centre - last_label >= LABEL_GAP:
                label, last_label = str(index + 1), centre
            count = int(self.exon_guides[index])
            coding = any(piece[2] for piece in pieces)
            exons.append({
                'number': index + 1,
                'blocks': blocks,
                'label': label,
                'label_x': round(centre, 2),
                'selected': selected_exon == index + 1,
                'title': 'Exon %d · %d bp%s · %d guide%s' % (
                    index + 1, high - low + 1, '' if coding else ' (untranslated)',
                    count, '' if count == 1 else 's'),
            })

        bins = []
        for index, sides in enumerate(self._bins.tolist()):
            both = [up + down for up, down in zip(*sides)]
            total = sum(both)
            if not total:
                continue
            x = MARGIN + index * bin_width
            segments = []
            for side, counts in enumerate(sides):
                edge = MIDLINE - 13 if side == 0 else MIDLINE + 13
                for name, value in zip(CLASSES, counts):
                    if not value:
                        continue
                    height = value / self.peak * BAR_MAX
                    y = edge - height if side == 0 else edge
                    segments.append({'x': round(x + 0.3, 2),
                                     'y': round(y, 2),
                                     'width': round(max(0.6, bin_width - 0.6), 2),
                                     'height': round(height, 2),
                                     'cls': name})
                    edge = edge - height if side == 0 else edge + height
            bins.append({
                'x': round(x, 2),
                'width': round(bin_width, 2),
                'segments': segments,
                'title': '%d guide%s: %s' % (
                    total, '' if total == 1 else 's',
                    ', '.join('%d %s' % (n, CLASS_LABELS[c])
                              for c, n in zip(CLASSES, both) if n)),
            })

        highlight = []
        if selected_substitution in self._substitution_rows:
            rows = self._substitution_rows[selected_substitution]
            for index in np.unique(self._row_bins[rows]):
                highlight.append({'x': round(MARGIN + index * bin_width, 2),
                                  'width': round(max(1.4, bin_width - 0.6), 2)})
        counts = [{'cls': name, 'label': CLASS_LABELS[name],
                   'count': self.class_counts.get(name, 0)}
                  for name in CLASSES if self.class_counts.get(name)]
        return {'exons': exons, 'bins': bins, 'highlight': highlight,
                'width': WIDTH, 'height': HEIGHT,
                'midline': MIDLINE, 'margin': MARGIN, 'counts': counts,
                'total': self.total}

    def matrix_view(self, selected=''):
        """The substitution matrix: 441 squares, whatever the gene."""
        # The scale spans the squares it shades: silent ones are drawn grey, so a
        # large synonymous count does not wash out the rest.
        shaded = [count for key, count in self.matrix_counts.items()
                  if count and key.split('-')[0] != key.split('-')[1]]
        low, high = min(shaded + [1]), max(shaded + [1])
        cells = []
        for down, source in enumerate(RESIDUES):
            for across, target in enumerate(RESIDUES):
                key = substitution_key(source, target)
                count = self.matrix_counts.get(key, 0)
                fill = None
                if count:
                    if source == target:
                        cls = 'sil'
                    else:
                        cls = 'fill'
                        fill = ramp_colour((count - low) / float(high - low)
                                           if high > low else 1.0)
                    sites = self.matrix_sites[key]
                    pathogenic = self.matrix_pathogenic.get(key, 0)
                    title = '%s to %s · %d guide%s at %d site%s%s' % (
                        ONE_LETTER[source], ONE_LETTER[target], count,
                        '' if count == 1 else 's', sites, '' if sites == 1 else 's',
                        ' · %d recreate a pathogenic variant' % pathogenic
                        if pathogenic else '')
                else:
                    cls, title = 'off', '%s to %s · not reachable with this editor' % (
                        ONE_LETTER[source], ONE_LETTER[target])
                cells.append({
                    'key': key,
                    'x': MATRIX_LEFT + across * CELL,
                    'y': MATRIX_TOP + down * CELL,
                    'count': count,
                    'cls': cls,
                    'fill': fill,
                    'selected': selected == key,
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
            'ramp': [{'offset': round(i / float(len(RAMP) - 1), 3), 'colour': colour}
                     for i, colour in enumerate(RAMP)],
            'low': low, 'high': high, 'shaded': bool(shaded),
        }

    def selection(self, key):
        """A one-line summary of the chosen square, for the caption."""
        source, target = key.split('-')
        return {
            'letters': '%s → %s' % (ONE_LETTER[source], ONE_LETTER[target]),
            'guides': self.matrix_counts[key],
            'sites': self.matrix_sites[key],
            'pathogenic': self.matrix_pathogenic.get(key, 0),
            'silent': source == target,
        }
