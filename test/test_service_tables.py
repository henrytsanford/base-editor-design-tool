"""Building a result table with its coverage panels.

Both front doors load a stored result through `load_table`, so the guards live here
rather than once per caller. The empty-exon case is the one that matters: it is the
difference between a table with no drawings above it and an IndexError out of
Geometry, and it is reachable from any transcript row the bundle cannot place.
"""
from bedesign import DESIGN_COLUMNS, tsv

from service.cachekey import result_key
from service.storage import LocalStorage
from service.tables import load_table

KEY = 'a' * 64
ROW = ['GGAGTTGGAAGCGCGTTACC', 'CCCC', 'MAP2K1', 'ENSG1', 'ENST1', '1', 'GRCh38']


class Shaped(object):
    """A reference that answers with whatever geometry a test wants."""

    def __init__(self, shape):
        self.shape = shape

    def geometry(self, transcript_id):
        return self.shape


def stored(tmp_path):
    storage = LocalStorage(str(tmp_path))
    padded = ROW + [''] * (len(DESIGN_COLUMNS) - len(ROW))
    storage.put(result_key(KEY, 'designs'), tsv.tsv_gz(DESIGN_COLUMNS, [padded]))
    return storage


def test_a_transcript_with_exons_gets_its_panels(tmp_path):
    references = Shaped({'exons': [[100, 200]], 'cds': [[120, 180]], 'strand': 1})
    table = load_table(stored(tmp_path), references, KEY, 'ENST1', 30)
    assert table.total == 1
    assert table.coverage is not None


def test_a_transcript_the_bundle_cannot_place_gets_a_table_and_no_panels(tmp_path):
    """None from the bundle is a drawable-nothing, not a failure."""
    table = load_table(stored(tmp_path), Shaped(None), KEY, 'ENST1', 30)
    assert table.total == 1
    assert table.coverage is None


def test_a_transcript_row_with_no_exons_gets_a_table_and_no_panels(tmp_path):
    """The guard that matters: Geometry indexes its first exon to lay out the map,
    so an empty list raises rather than drawing nothing."""
    references = Shaped({'exons': [], 'cds': [], 'strand': 1})
    table = load_table(stored(tmp_path), references, KEY, 'ENST1', 30)
    assert table.total == 1
    assert table.coverage is None


def test_a_run_with_no_transcript_is_never_given_geometry(tmp_path):
    """A pasted sequence has no transcript to ask the bundle about."""
    def explode(transcript_id):
        raise AssertionError('the bundle should not be consulted')

    references = Shaped(None)
    references.geometry = explode
    table = load_table(stored(tmp_path), references, KEY, '', 30)
    assert table.coverage is None
