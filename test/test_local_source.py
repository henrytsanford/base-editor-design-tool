"""Differential test: the local bundle must answer exactly as Ensembl does.

This is the test the whole local-mirror design rests on. Everything else checks
that the code is self-consistent; only this checks that the reimplementation of
Ensembl's coordinate conventions is *right* — stop codons folded into the CDS,
blocks in transcript order, flanks clamped at contig edges, minus-strand
sequence reverse-complemented.

Opt in, because it queries the live REST API and needs a built bundle:

    ENSEMBL_TESTS=1 pytest test/test_local_source.py
"""
import pytest

import transcript_source as ts

pytestmark = [pytest.mark.ensembl, pytest.mark.bundle]

# One transcript per way the reimplementation could go wrong. Names are from
# Ensembl release 116; the IDs are stable across releases.
PANEL = [
    ("ENST00000294952", "PPP1R21-202, + strand, 22 exons (the original bug report)"),
    ("ENST00000369550", "DKC1-201, + strand, chrX, 15 exons"),
    ("ENST00000269305", "TP53-201, - strand, 11 exons"),
    ("ENST00000352993", "BRCA1-201, - strand, 22 exons"),
    ("ENST00000325404", "SOX2-201, + strand, single exon"),
    ("ENST00000371222", "JUN-201, - strand, single exon"),
    ("ENST00000508832", "MALAT1-201, non-coding: no CDS mapping at all"),
    ("ENST00000320065", "OR2G2-201, + strand, no UTR: CDS fills the transcript"),
    ("ENST00000334857", "OR10J5-201, - strand, no UTR"),
]

# Every field get_tr_info and get_utrs read out of the lookup response.
LOOKUP_FIELDS = ("display_name", "assembly_name", "strand", "seq_region_name",
                 "Parent", "start", "end")


@pytest.fixture(scope="module")
def sources(refdata):
    return ts.LocalSource(ts.find_bundle(refdata)), ts.EnsemblRestSource()


@pytest.fixture(scope="module")
def cache():
    return {}


def both(sources, cache, tr):
    """Answers from both sources, fetched once per transcript."""
    if tr not in cache:
        local, rest = sources
        span = rest.lookup(tr)
        # get_tr_info computes the /map/cds length this way, for either strand
        length = span["end"] - span["start"]
        cache[tr] = {
            "rest": {
                "lookup": span,
                "cds_mappings": rest.cds_mappings(tr, length),
                "genomic": rest.genomic_sequence(tr),
                "protein": rest.protein_sequence(tr),
                "cds": rest.cds_sequence(tr),
            },
            "local": {
                "lookup": local.lookup(tr),
                "cds_mappings": local.cds_mappings(tr, length),
                "genomic": local.genomic_sequence(tr),
                "protein": local.protein_sequence(tr),
                "cds": local.cds_sequence(tr),
            },
        }
    return cache[tr]["local"], cache[tr]["rest"]


@pytest.mark.parametrize("tr,description", PANEL, ids=[t for t, _ in PANEL])
def test_lookup_agrees(sources, cache, tr, description):
    local, rest = both(sources, cache, tr)
    for field in LOOKUP_FIELDS:
        assert local["lookup"][field] == rest["lookup"][field], \
            f"{field} differs for {description}"


@pytest.mark.parametrize("tr,description", PANEL, ids=[t for t, _ in PANEL])
def test_exons_agree(sources, cache, tr, description):
    """Order matters: get_utrs reads Exon[0] as the transcript's first exon."""
    local, rest = both(sources, cache, tr)
    assert [(e["start"], e["end"]) for e in local["lookup"]["Exon"]] == \
           [(e["start"], e["end"]) for e in rest["lookup"]["Exon"]], description


@pytest.mark.parametrize("tr,description", PANEL, ids=[t for t, _ in PANEL])
def test_cds_mappings_agree(sources, cache, tr, description):
    """Catches a mis-handled stop codon, which would shift every coordinate 3 bp."""
    local, rest = both(sources, cache, tr)
    if rest["cds_mappings"] is None:
        assert local["cds_mappings"] is None, description
        return
    assert [(m["start"], m["end"], m["strand"]) for m in local["cds_mappings"]] == \
           [(m["start"], m["end"], m["strand"]) for m in rest["cds_mappings"]], description


@pytest.mark.parametrize("tr,description", PANEL, ids=[t for t, _ in PANEL])
def test_genomic_sequence_agrees(sources, cache, tr, description):
    """Covers strand orientation, flank placement and contig clamping at once."""
    local, rest = both(sources, cache, tr)
    assert local["genomic"] == rest["genomic"], description


@pytest.mark.parametrize("tr,description", PANEL, ids=[t for t, _ in PANEL])
def test_protein_and_cds_sequences_agree(sources, cache, tr, description):
    local, rest = both(sources, cache, tr)
    assert local["protein"] == rest["protein"], description
    assert local["cds"] == rest["cds"], description


def test_unknown_transcript_raises_from_both_sources(sources):
    local, rest = sources
    with pytest.raises(ts.TranscriptNotFound):
        local.lookup("ENST00000000000")
    with pytest.raises(ts.TranscriptNotFound):
        rest.lookup("ENST00000000000")
