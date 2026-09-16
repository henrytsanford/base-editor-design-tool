"""Unit tests for the reference-data layer.

LocalSource is exercised against a small synthetic bundle built in a tmpdir, so
these run offline and pin down the semantics that have to match the REST API
exactly: transcript-order blocks, the stop codon folded into the CDS, flanks
clamped at contig edges, and minus-strand sequence reverse-complemented.
"""
import gzip
import json
import os
import shutil
import sqlite3

import pytest
from Bio import bgzf

import bedesign.transcript_source as ts
from bedesign import engine
from bedesign.transcript_source import TRANSCRIPT_COLUMNS
from build_reference import (SCHEMA, merge_stop_codon,
                             order_blocks, parse_attributes, read_gtf)


# --------------------------------------------------------------- GTF handling

def test_parse_attributes():
    attrs = parse_attributes(
        'gene_id "ENSG1"; transcript_id "ENST1"; exon_number "3"; '
        'transcript_name "ABC-201"; tag "basic";\n')
    assert attrs['gene_id'] == 'ENSG1'
    assert attrs['exon_number'] == '3'
    assert attrs['transcript_name'] == 'ABC-201'


def test_parse_attributes_collects_repeated_tags():
    """A GTF line repeats `tag`; a dict keyed by name would keep only the last."""
    attrs = parse_attributes(
        'gene_id "ENSG1"; transcript_id "ENST1"; tag "CCDS"; ccds_id "CCDS152"; '
        'tag "gencode_basic"; tag "MANE_Select"; tag "Ensembl_canonical";\n')
    assert attrs['tag'] == ['CCDS', 'gencode_basic', 'MANE_Select',
                            'Ensembl_canonical']
    # single-valued keys are untouched, including one sitting between two tags
    assert attrs['ccds_id'] == 'CCDS152'
    assert attrs['gene_id'] == 'ENSG1'


def test_parse_attributes_without_tags_has_no_tag_key():
    attrs = parse_attributes('gene_id "ENSG1"; transcript_id "ENST1";\n')
    assert attrs.get('tag', []) == []


GTF_ATTRS = ('gene_id "ENSG1"; transcript_id "%s"; gene_name "G"; '
             'transcript_name "G-20%s"; transcript_biotype "protein_coding"; %s')


def write_gtf(path, lines):
    with gzip.open(str(path), 'wt') as fh:
        fh.write('#!genome-version GRCh38\n')
        for attrs in lines:
            fh.write('\t'.join(['7', 'ensembl', 'transcript', '101', '300',
                                '.', '+', '.', attrs]) + '\n')


def test_read_gtf_flags_mane_and_canonical(tmp_path):
    gtf = tmp_path / 'test.gtf.gz'
    write_gtf(gtf, [
        GTF_ATTRS % ('ENSTMANE', '1',
                     'tag "CCDS"; tag "MANE_Select"; tag "Ensembl_canonical";'),
        GTF_ATTRS % ('ENSTCANON', '2', 'tag "basic"; tag "Ensembl_canonical";'),
        GTF_ATTRS % ('ENSTPLAIN', '3', 'tag "basic";'),
    ])
    transcripts, assembly = read_gtf(str(gtf))
    assert assembly == 'GRCh38'
    assert (transcripts['ENSTMANE']['mane_select'],
            transcripts['ENSTMANE']['ensembl_canonical']) == (1, 1)
    assert (transcripts['ENSTCANON']['mane_select'],
            transcripts['ENSTCANON']['ensembl_canonical']) == (0, 1)
    assert (transcripts['ENSTPLAIN']['mane_select'],
            transcripts['ENSTPLAIN']['ensembl_canonical']) == (0, 0)


def test_order_blocks_is_transcript_order():
    blocks = [(2, 200, 250), (1, 100, 150), (3, 300, 350)]
    assert order_blocks(blocks, 1) == [(100, 150), (200, 250), (300, 350)]
    # minus strand runs 5'->3' down the chromosome
    assert order_blocks(blocks, -1) == [(300, 350), (200, 250), (100, 150)]


def test_stop_codon_extends_the_last_cds_block():
    # Ensembl's CDS coordinates include the stop codon; the GTF stores it apart.
    assert merge_stop_codon([(100, 150), (200, 250)], [(251, 253)], 1) \
        == [(100, 150), (200, 253)]
    assert merge_stop_codon([(300, 350), (100, 150)], [(97, 99)], -1) \
        == [(300, 350), (97, 150)]


def test_split_stop_codon_becomes_its_own_block():
    """A stop codon interrupted by an intron does not abut the last CDS block."""
    assert merge_stop_codon([(100, 150)], [(400, 402)], 1) \
        == [(100, 150), (400, 402)]


def test_stop_codon_absent_leaves_cds_alone():
    assert merge_stop_codon([(100, 150)], [], 1) == [(100, 150)]


# ------------------------------------------------------------ a tiny bundle

GENOME = ''.join(['A' * 100, 'C' * 100, 'G' * 100, 'T' * 100])  # 400 bp, contig "7"

# plus strand, two exons, CDS 121-150 + 201-260 (90 bp incl. stop codon)
PLUS = {
    'transcript_id': 'ENST0000PLUS', 'display_name': 'PLUS-201',
    'gene_id': 'ENSG0000PLUS', 'gene_name': 'PLUS', 'biotype': 'protein_coding',
    'seq_region': '7', 'strand': 1, 'start': 101, 'end': 300,
    'exons': [[101, 150], [201, 300]], 'cds': [[121, 150], [201, 260]],
    'cds_sequence': 'ATG' + 'AAA' * 28 + 'TGA', 'protein_sequence': 'M' + 'K' * 28,
    'mane_select': 1, 'ensembl_canonical': 1,
}
# minus strand, same span. Canonical for its gene but not MANE, the common case
# for genes RefSeq and Ensembl have not agreed on.
MINUS = dict(PLUS, transcript_id='ENST0000MINUS', display_name='MINUS-201',
             gene_id='ENSG0000MINUS', gene_name='MINUS', strand=-1,
             exons=[[201, 300], [101, 150]], cds=[[201, 260], [121, 150]],
             mane_select=0, ensembl_canonical=1)
# runs to the very start of the contig, so a 40 bp flank cannot fit
EDGE = dict(PLUS, transcript_id='ENST0000EDGE', display_name='EDGE-201',
            gene_id='ENSG0000EDGE', gene_name='EDGE', start=10, end=60,
            exons=[[10, 60]], cds=[[10, 60]])
# non-coding: no CDS at all
NONCODING = dict(PLUS, transcript_id='ENST0000NC', display_name='NC-201',
                 gene_id='ENSG0000NC', gene_name='NC', biotype='lncRNA',
                 cds=None, cds_sequence=None, protein_sequence=None)

RECORDS = [PLUS, MINUS, EDGE, NONCODING]


@pytest.fixture(scope='module')
def bundle(tmp_path_factory):
    path = tmp_path_factory.mktemp('refdata') / 'ensembl-999'
    path.mkdir()
    db = sqlite3.connect(str(path / 'transcripts.db'))
    db.executescript(SCHEMA)
    insert = ('INSERT INTO transcript (%s) VALUES (%s)'
              % (','.join(TRANSCRIPT_COLUMNS), ','.join('?' * len(TRANSCRIPT_COLUMNS))))
    for r in RECORDS:
        db.execute(insert, (
            r['transcript_id'], r['display_name'], r['gene_id'], r['gene_name'],
            r['biotype'], r['seq_region'], r['strand'], r['start'], r['end'],
            json.dumps(r['exons']),
            json.dumps(r['cds']) if r['cds'] else None,
            r['cds_sequence'], r['protein_sequence'],
            r['mane_select'], r['ensembl_canonical']))
    db.executemany('INSERT INTO meta VALUES (?,?)',
                   [('release', '999'), ('assembly', 'GRCh38')])
    db.commit()
    db.close()
    fa = str(path / 'genome.fa.bgz')
    with bgzf.BgzfWriter(fa, 'wb') as w:
        w.write(('>7 dna:chromosome\n'
                 + '\n'.join(GENOME[i:i + 60] for i in range(0, len(GENOME), 60))
                 + '\n').encode())
    from pyfaidx import Fasta
    Fasta(fa)
    return str(path)


@pytest.fixture
def local(bundle):
    return ts.LocalSource(bundle)


def test_lookup_shape_matches_the_rest_response(local):
    info = local.lookup('ENST0000PLUS')
    # exactly the fields get_tr_info and get_utrs read
    assert info['display_name'] == 'PLUS-201'
    assert info['assembly_name'] == 'GRCh38'
    assert info['strand'] == 1
    assert info['seq_region_name'] == '7'
    assert info['Parent'] == 'ENSG0000PLUS'
    assert info['start'] == 101 and info['end'] == 300
    assert info['Exon'] == [{'start': 101, 'end': 150}, {'start': 201, 'end': 300}]


def test_lookup_reports_mane_and_canonical(local):
    """What the gene search preselects. REST has no equivalent field."""
    plus = local.lookup('ENST0000PLUS')
    assert plus['gene_name'] == 'PLUS'
    assert plus['mane_select'] is True
    assert plus['ensembl_canonical'] is True
    # canonical for its gene, but not the MANE Select transcript
    minus = local.lookup('ENST0000MINUS')
    assert minus['mane_select'] is False
    assert minus['ensembl_canonical'] is True


def test_a_bundle_without_the_tag_columns_is_refused(bundle, tmp_path):
    """A bundle is a build artifact, so an old one is rebuilt, not worked around.

    Answering lookups with mane_select=False for every transcript would make the
    gene search silently preselect nothing, which reads as a bug in the caller.
    """
    old = str(tmp_path / 'ensembl-998')
    shutil.copytree(bundle, old)
    db = sqlite3.connect(os.path.join(old, ts.TRANSCRIPTS_DB))
    kept = [c for c in TRANSCRIPT_COLUMNS
            if c not in ('mane_select', 'ensembl_canonical')]
    db.executescript(
        'CREATE TABLE t AS SELECT %s FROM transcript;'
        'DROP TABLE transcript;'
        'ALTER TABLE t RENAME TO transcript;' % ','.join(kept))
    db.commit()
    db.close()
    with pytest.raises(ts.BundleNotFound, match='rebuilt'):
        ts.LocalSource(old)


def test_lookup_exons_are_in_transcript_order(local):
    """get_utrs treats Exon[0] as the transcript's first exon, not the leftmost."""
    assert local.lookup('ENST0000MINUS')['Exon'][0] == {'start': 201, 'end': 300}


def test_unknown_transcript_raises_instead_of_exiting(local):
    with pytest.raises(ts.TranscriptNotFound) as excinfo:
        local.lookup('ENST0000NOPE')
    assert 'ENST0000NOPE' in str(excinfo.value)


def test_cds_mappings_are_transcript_ordered_blocks(local):
    m = local.cds_mappings('ENST0000PLUS', 199)
    assert [(b['start'], b['end']) for b in m] == [(121, 150), (201, 260)]


def test_cds_mappings_none_when_transcript_has_no_cds(local):
    assert local.cds_mappings('ENST0000NC', 199) is None


def test_cds_mappings_none_rather_than_empty(local):
    """get_exons indexes mappings[0], so an empty list must not reach it."""
    assert local.cds_mappings('ENST0000PLUS', 0) is None


def test_cds_mappings_honour_the_requested_length(local):
    """Callers pass the genomic span, which for a single-exon CDS can fall one
    base short of the CDS; the REST API stops at the requested length."""
    m = local.cds_mappings('ENST0000PLUS', 40)
    assert [(b['start'], b['end']) for b in m] == [(121, 150), (201, 210)]


def test_cds_mappings_trim_from_the_three_prime_end_on_minus_strand(local):
    m = local.cds_mappings('ENST0000MINUS', 70)
    assert [(b['start'], b['end']) for b in m] == [(201, 260), (141, 150)]


def test_genomic_sequence_adds_flanks(local):
    seq = local.genomic_sequence('ENST0000PLUS', flank=40)
    assert seq == GENOME[60:340]
    assert len(seq) == (300 - 101 + 1) + 80


def test_genomic_sequence_is_reverse_complemented_on_minus_strand(local):
    plus = local.genomic_sequence('ENST0000PLUS', flank=40)
    minus = local.genomic_sequence('ENST0000MINUS', flank=40)
    assert minus == ts.reverse_complement(plus)


def test_genomic_sequence_clamps_at_the_contig_edge(local):
    seq = local.genomic_sequence('ENST0000EDGE', flank=40)
    # the transcript starts at 10, so only 9 bp of 5' flank exist
    assert seq == GENOME[0:100]


def test_cds_and_protein_sequences(local):
    assert local.cds_sequence('ENST0000PLUS').startswith('ATG')
    assert local.protein_sequence('ENST0000PLUS').startswith('M')


def test_missing_sequences_come_back_empty_not_none(local):
    assert local.cds_sequence('ENST0000NC') == ''
    assert local.protein_sequence('ENST0000NC') == ''


# ------------------------------------------------------------ bundle discovery

def test_find_bundle_accepts_the_bundle_itself(bundle):
    assert ts.find_bundle(bundle) == bundle


def test_find_bundle_picks_the_highest_release(bundle):
    parent = os.path.dirname(bundle)
    older = os.path.join(parent, 'ensembl-9')
    os.makedirs(older, exist_ok=True)
    for name in ts.BUNDLE_FILES:
        open(os.path.join(older, name), 'a').close()
    assert ts.find_bundle(parent) == bundle


def test_missing_bundle_is_reported_not_traced(tmp_path):
    with pytest.raises(ts.BundleNotFound):
        ts.find_bundle(str(tmp_path))
    with pytest.raises(ts.BundleNotFound):
        ts.LocalSource(str(tmp_path))


def test_a_half_built_bundle_is_not_usable(tmp_path):
    """build_reference.py --stage db leaves a directory with no genome in it."""
    partial = tmp_path / 'ensembl-1'
    partial.mkdir()
    (partial / 'transcripts.db').touch()
    with pytest.raises(ts.BundleNotFound) as excinfo:
        ts.LocalSource(str(partial))
    assert 'genome.fa.bgz' in str(excinfo.value)
    with pytest.raises(ts.BundleNotFound):
        ts.find_bundle(str(tmp_path))


# ------------------------------------------------- the engine takes a source

class StubSource:
    """Anything with these methods will do; the engine never names a class."""

    def genomic_sequence(self, tr, flank=40):
        return 'ACGT'

    def protein_sequence(self, tr):
        return 'MK'

    def cds_sequence(self, tr):
        return 'ATGAAA'


def test_getters_read_the_source_they_are_given():
    """The source is an argument, not module state, so a caller can hold several.

    That is what lets the web service run designs in a process pool against one
    bundle and one ClinVar database per worker.
    """
    source = StubSource()
    assert engine.get_tr_sequence(source, 'ENST1') == 'ACGT'
    assert engine.get_pro_sequence(source, 'ENST1') == 'MK'
    assert engine.get_cds_sequence(source, 'ENST1') == 'ATGAAA'


def test_get_source_rejects_an_unknown_name():
    with pytest.raises(ValueError):
        ts.get_source('sqlite3')
