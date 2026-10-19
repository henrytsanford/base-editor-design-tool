"""Ambiguity codes in a transcript's genomic span.

GRCh38 hard-masks occasional bases as 'N'. They sit in the genomic sequence the
design loop scans, not in the CDS, so they cost the guides whose context window
covers them -- and nothing else. CD4 (ENST00000011653) carries two of them in an
intron, thousands of bases from any exon, and used to produce no designs at all.
"""
import pytest

from bedesign.engine import (DesignParams, check_sequences, count_ambiguous,
                             design_transcript)
from bedesign.transcript_source import get_source

# Two 'N' bases in intron 3, both far outside the 30 bp the design loop reaches.
CD4 = "ENST00000011653"
# Single exon, 450 bp of genomic span: small enough to design twice in a test.
SMALL = "ENST00000416044"  # KRTAP25-1


class MaskedSource(object):
    """A source whose genomic sequence has one base replaced by 'N'.

    Wraps rather than subclasses, so it stands in for any TranscriptSource and
    changes only the sequence the design loop scans -- the CDS, protein and
    coordinates come through untouched, which is the shape of a real hard mask.
    """

    def __init__(self, inner, index):
        self.inner = inner
        self.index = index

    def __getattr__(self, name):
        return getattr(self.inner, name)

    def genomic_sequence(self, tr, flank=40):
        seq = self.inner.genomic_sequence(tr, flank)
        return seq[:self.index] + "N" + seq[self.index + 1:]


@pytest.fixture(scope="module")
def local_source(refdata):
    return get_source("local", refdata)


def test_check_sequences_rejects_ambiguity_by_default():
    assert check_sequences("ACTGN") != ""
    assert check_sequences("ACTG") == ""


def test_check_sequences_can_allow_ambiguity():
    assert check_sequences("ACTGNRYKM", allow_ambiguous=True) == ""
    # Still not a nucleotide sequence, whatever the caller tolerates
    assert check_sequences("ACTGX", allow_ambiguous=True) != ""
    assert check_sequences("MKRAPLE", allow_ambiguous=True) != ""


def test_count_ambiguous():
    assert count_ambiguous("ACTGACTG") == 0
    assert count_ambiguous("ACTNGnACTG") == 2


@pytest.mark.bundle
def test_intronic_n_still_designs(local_source):
    """CD4's masked intronic bases must not cost the whole transcript."""
    designs, errors, _ = design_transcript(
        local_source, None, CD4, DesignParams.from_preset("ABE8e-SpRY"))
    assert designs
    assert errors == []


@pytest.mark.bundle
def test_only_overlapping_guides_are_lost(local_source):
    """An 'N' inside the designed region costs the guides that span it.

    Comparing a masked run against an unmasked one pins both halves of the fix:
    every surviving row is unchanged, and the rows that go missing are only the
    handful whose context window could reach the masked base.
    """
    params = DesignParams(pam="NGG", window="4-8", sg_len=20, edit="C-T")
    index = len(local_source.genomic_sequence(SMALL)) // 2

    designs, errors, _ = design_transcript(local_source, None, SMALL, params)
    masked, masked_errors, _ = design_transcript(
        MaskedSource(local_source, index), None, SMALL, params)

    assert errors == [] and masked_errors == []
    dropped = [row for row in designs if row not in masked]
    assert len(masked) == len(designs) - len(dropped)
    assert dropped, "the masked base was outside the designed region"
    # A base reaches at most one context window per start position, on either strand
    context_len = params.sg_len + len(params.pam) + 8
    assert len(dropped) <= 2 * context_len
