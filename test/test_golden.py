"""Golden output fixtures: the engine's output, frozen.

Runs the CLI over a panel of transcripts that between them reach every
structural branch in the design code, and compares all three output files byte
for byte against committed fixtures. Any change to any output row fails here.

The comparison is whole-file, with no dropped columns -- unlike
assert_output_files_equal in test_design.py, which excludes three columns to
match an older reference dataset. `Clinical significance` is covered here.

To rewrite the fixtures after an intentional output change:

    pytest test/test_golden.py --regenerate-golden

and read the diff before committing it.
"""
import gzip
import os

import pytest

from bedesign import ANNOTATIONS_FILE, DESIGNS_FILE, ERRORS_FILE

GOLDEN_DIR = 'test/golden'

# The files a run writes, as (fixture name, filename template). Taken from the
# engine so a rename shows up here as a rename, not as a missing fixture.
DESIGN_OUTPUTS = (('sgrna_designs', DESIGNS_FILE), ('error_report', ERRORS_FILE))
ANNOTATION_OUTPUT = ('clinvar_annotations', ANNOTATIONS_FILE)
# Nucleotide input has no gene to look up, so it writes no annotations file.
TID_OUTPUTS = DESIGN_OUTPUTS + (ANNOTATION_OUTPUT,)
NUC_OUTPUTS = DESIGN_OUTPUTS

# Transcripts chosen to reach every branch in design_sgrnas and
# get_context_for_trans, kept small so the fixtures stay well under a megabyte.
PANEL = [
    ('map2k1', 'ENST00000307102', 'MAP2K1'),    # plus strand, 11 exons, 659 ClinVar variants
    ('isy1', 'ENST00000393295', 'ISY1'),        # minus strand, 11 exons, a few variants
    ('krtap25_1', 'ENST00000416044', 'KRTAP25-1'),  # single exon with UTR, minus strand
    ('hgs', 'ENST00000678176', 'HGS'),          # no UTR, 2 exons, CDS not a multiple of 3
    ('malat1', 'ENST00000620902', 'MALAT1'),    # non-coding: no CDS at all
]


# Run names are prefixed so a run can never collide with -- and so the fixture
# teardown can never delete -- an output directory someone made by hand.
RUN_PREFIX = 'goldenrun_'


def golden_path(name, output):
    return os.path.join(GOLDEN_DIR, name, output + '.txt.gz')


def read_golden(name, output):
    with gzip.open(golden_path(name, output), 'rt') as fh:
        return fh.read()


def write_golden(name, output, text):
    path = golden_path(name, output)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # mtime=0 so regenerating unchanged output produces an identical file
    with gzip.GzipFile(path, 'wb', compresslevel=9, mtime=0) as fh:
        fh.write(text.encode())


def check(name, output_dir, output_name, outputs, regenerate):
    """Compares, or rewrites, every output the run was expected to write.

    Every file in `outputs` must exist: a run that silently stopped writing one
    should fail here, not be skipped over.
    """
    for output, template in outputs:
        filename = template % output_name if '%s' in template else template
        produced = os.path.join(output_dir, filename)
        assert os.path.exists(produced), 'run wrote no %s' % produced
        with open(produced) as fh:
            text = fh.read()
        if regenerate:
            write_golden(name, output, text)
            continue
        assert os.path.exists(golden_path(name, output)), (
            'no golden fixture for %s/%s; run pytest --regenerate-golden'
            % (name, output))
        assert text == read_golden(name, output), (
            '%s/%s differs from its golden fixture. If the change is intended, '
            'rerun with --regenerate-golden and review the diff.' % (name, output))


@pytest.fixture
def regenerate(request):
    return request.config.getoption('--regenerate-golden')


@pytest.mark.bundle
@pytest.mark.parametrize('name,transcript,gene', PANEL, ids=[p[0] for p in PANEL])
def test_transcript_matches_golden(name, transcript, gene,
                                   run_design, tmp_path, refdata, regenerate):
    """Every transcript in the panel reproduces its frozen output exactly."""
    path = tmp_path / (name + '.txt')
    path.write_text('Transcript ID\tGene Symbol\n%s\t%s\n' % (transcript, gene))
    run_name = RUN_PREFIX + name
    # --edit all runs both the C-T and A-G passes, which is what the service will do
    _, output_dir = run_design(str(path), 'tid', run_name,
                               '--source', 'local', '--refdata', refdata,
                               '--edit', 'all')
    check(name, output_dir, run_name, TID_OUTPUTS, regenerate)


def test_fasta_matches_golden(run_design, regenerate):
    """The nucleotide path, which needs neither a bundle nor the network."""
    _, output_dir = run_design('Sample_data/GFP.fasta', 'nuc', RUN_PREFIX + 'gfp',
                               '--edit', 'all')
    check('gfp', output_dir, RUN_PREFIX + 'gfp', NUC_OUTPUTS, regenerate)


@pytest.mark.bundle
def test_golden_designs_carry_clinical_significance():
    """ClinVar annotation reaches the designs file, not just the annotations file.

    MAP2K1 is the densest gene in the ClinVar fixture, so if its frozen output
    carries no significance at all, the annotation has silently broken.
    """
    designs = read_golden('map2k1', 'sgrna_designs').splitlines()
    header = designs[0].split('\t')
    column = header.index('Clinical significance')
    reported = {row.split('\t')[column].replace('None', '').strip('; ')
                for row in designs[1:]}
    assert any(reported - {''}), 'no guide in the MAP2K1 fixture carries a significance'
