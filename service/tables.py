"""Building a parsed result table, with its coverage panels attached.

Both front doors need the same thing from a stored result: the designs frame, plus
the two coverage panels when the bundle can describe the transcript's shape. That is
one job with one set of guards, so it lives here rather than once per caller -- the
copies drifted the first time there were two.
"""
from bedesign import DESIGN_COLUMNS

from . import coverage
from .cachekey import result_key
from .results import ResultTable


def load_table(storage, references, key, transcript_id, intron_buffer):
    """The designs table for a stored result, with coverage where it is drawable.

    A transcript the bundle cannot place -- no row, or a row carrying no exons --
    leaves `coverage` as None, and the caller shows a table with no panels above it.
    An empty exon list has to be caught here rather than passed on: Geometry indexes
    its first exon while laying out the drawing.
    """
    table = ResultTable(storage.get(result_key(key, 'designs')), DESIGN_COLUMNS)
    shape = references.geometry(transcript_id) if transcript_id else None
    if shape is not None and shape['exons']:
        table.coverage = coverage.Coverage(
            table, coverage.Geometry(shape['exons'], shape['cds'], shape['strand'],
                                     intron_buffer))
    return table
