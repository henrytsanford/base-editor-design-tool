import os
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest

import build_clinvar
from base_editing_guide_designs import get_aa_map

SAMPLE_DATA_DIR = "Sample_data"
SAMPLE_REFERENCE = f"{SAMPLE_DATA_DIR}/sample_19-05-07-14-35-42"

# Real ClinVar rows for three of the sample genes (ISY1, PSMB5, MAP2K1), in
# variant_summary.txt format.
CLINVAR_FIXTURE = "test/data/variant_summary_sample.txt.gz"
CLINVAR_GENE = ("ENST00000307102", "MAP2K1")  # 659 variants, densest of the three

# ENST00000334810 (ADGRD2) is in the sample input but was retired from Ensembl
# after the reference output was made. The REST API still resolves retired IDs
# from its archive; a release bundle contains only what is current in that
# release, so the local source cannot design for it.
RETIRED_IN_CURRENT_RELEASE = ["ENST00000334810"]


@pytest.fixture(scope="session")
def clinvar_db(tmp_path_factory):
    """A ClinVar database over the sample genes, built by the real builder.

    Building it here rather than checking in a .db keeps the fixture honest: a
    change to build_clinvar.py that broke the schema or the column mapping would
    break these tests too.
    """
    path = tmp_path_factory.mktemp("clinvar") / "clinvar-test.db"
    build_clinvar.build(CLINVAR_FIXTURE, str(path))
    return str(path)


@pytest.fixture
def run_design(clinvar_db):
    """Runs the script, returning (result, output folder); removes the output afterwards."""
    folders = []

    def run(input_file, input_type, output_name, *extra, check=True):
        cmd = [sys.executable, "base_editing_guide_designs.py",
               "--input-file", input_file, "--input-type", input_type,
               "--clinvar-db", clinvar_db, "--pam", "NGG", "--intron-buffer", "30",
               "--edit", "C-T", "--output-name", output_name, "--sg-len", "20", *extra]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if check:
            assert result.returncode == 0, (
                f"Script failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\n"
                f"STDERR: {result.stderr}"
            )
        pattern = re.compile(re.escape(output_name) + r"_\d\d(-\d\d){5}$")
        found = [d for d in os.listdir(".") if pattern.match(d)]
        folders.extend(found)
        return result, found[0] if found else None

    yield run
    for d in folders:
        shutil.rmtree(d, ignore_errors=True)


def assert_output_files_equal(output_dir, sample_data_results_dir, output_name,
                              exclude=(), check_readme=True):
    """Compares a run against the checked-in reference output, ignoring the
    Ensembl transcript IDs in `exclude`."""
    output_file = f"sgrna_designs_{output_name}"
    with open(f"{output_dir}/{output_file}.txt") as f1, open(
        f"{sample_data_results_dir}/{output_file}.txt"
    ) as f2:
        # don't test columns with known differences:
        # Ensembl Gene ID missing from test file
        # Clinical significance depends on the ClinVar release the run used
        # BsmBI flag logic changed from test file
        df1 = pd.read_csv(f1, sep="\t").drop(
            columns=["Ensembl Gene ID", "Clinical significance", "BsmBI flag"])
        df2 = pd.read_csv(f2, sep="\t").drop(
            columns=["Clinical significance", "BsmBI flag"])
        if exclude:
            column = "Ensembl transcript ID"
            df1 = df1[~df1[column].isin(exclude)].reset_index(drop=True)
            df2 = df2[~df2[column].isin(exclude)].reset_index(drop=True)
        f1_contents = df1.to_csv(sep="\t", index=False)
        f2_contents = df2.to_csv(sep="\t", index=False)

        # test file uses single character amino acid codes
        aa_map = get_aa_map()
        aa_map["Ter"] = "*"
        for key, value in aa_map.items():
            f1_contents = f1_contents.replace(key, value)

        assert f1_contents == f2_contents
    if not check_readme:
        return
    with open(f"{output_dir}/README.txt") as f1, open(
        f"{sample_data_results_dir}/README.txt"
    ) as f2:
        # final two lines contain timestamps, don't check them
        assert f1.readlines()[:5] == f2.readlines()[:5]


def test_fasta_input(run_design):
    _, output_dir = run_design(f"{SAMPLE_DATA_DIR}/GFP.fasta", "nuc", "GFP")
    assert_output_files_equal(output_dir, f"{SAMPLE_DATA_DIR}/GFP_19-05-07-14-41-15", "GFP")


@pytest.mark.ensembl
def test_tid_input(run_design):
    _, output_dir = run_design(f"{SAMPLE_DATA_DIR}/sample_input.txt", "tid", "sample")
    assert_output_files_equal(output_dir, SAMPLE_REFERENCE, "sample")


@pytest.fixture
def current_sample_input(tmp_path):
    """The sample input minus transcripts retired since it was written."""
    path = tmp_path / "sample_input_current.txt"
    with open(f"{SAMPLE_DATA_DIR}/sample_input.txt") as fh:
        lines = fh.readlines()
    keep = [lines[0]] + [l for l in lines[1:]
                         if l.split("\t")[0].strip() not in RETIRED_IN_CURRENT_RELEASE]
    path.write_text("".join(keep))
    return str(path)


@pytest.mark.bundle
def test_tid_input_local_source(run_design, current_sample_input, refdata):
    """The local mirror must reproduce the reference output, with no network.

    If it does, the mirror is faithful for everything the designs depend on: exon
    and CDS coordinates, strand, flanks, and the CDS and protein sequences.
    """
    _, output_dir = run_design(current_sample_input, "tid", "sample",
                               "--source", "local", "--refdata", refdata)
    # the README names the temporary input file, so it cannot match the reference
    assert_output_files_equal(output_dir, SAMPLE_REFERENCE, "sample",
                              exclude=RETIRED_IN_CURRENT_RELEASE, check_readme=False)


@pytest.mark.bundle
def test_clinvar_annotates_a_known_gene(run_design, tmp_path, refdata):
    """A gene with hundreds of ClinVar SNVs must come back with annotations.

    Empty annotation is invisible to every other test here: get_snps writes a row
    per edit whether or not a SNP matched, and the golden comparison drops the
    Clinical significance column. So assert on matches, not on row counts.
    """
    transcript, gene = CLINVAR_GENE
    path = tmp_path / "clinvar_gene.txt"
    path.write_text(f"Transcript ID\tGene Symbol\n{transcript}\t{gene}\n")
    _, output_dir = run_design(str(path), "tid", "clinvar_gene",
                               "--source", "local", "--refdata", refdata)

    annotations = pd.read_csv(f"{output_dir}/clinvar_annotations_clinvar_gene.txt",
                              sep="\t")
    matched = annotations[annotations["SNP name"].notna()]
    assert len(matched) > 0, f"no ClinVar SNP matched any edit in {gene}"
    assert matched["SNP clinical significance"].notna().all()

    designs = pd.read_csv(f"{output_dir}/sgrna_designs_clinvar_gene.txt", sep="\t")
    # get_clinical_sig joins 'None' per unmatched edit, so a working run is one
    # where at least some guides carry something other than 'None'.
    reported = designs["Clinical significance"].fillna("").str.replace("None", "")
    assert (reported.str.strip(";").str.strip() != "").any()


@pytest.mark.bundle
def test_retired_transcript_is_reported_clearly(run_design, tmp_path, refdata):
    """A bundle holds one Ensembl release, so IDs retired since then are absent.

    Users paste old IDs from papers and spreadsheets, so the message has to say
    what to do rather than just fail.
    """
    path = tmp_path / "retired.txt"
    path.write_text("Transcript ID\tGene Symbol\nENST00000334810\tADGRD2\n")
    result, _ = run_design(str(path), "tid", "sample_retired",
                           "--source", "local", "--refdata", refdata, check=False)
    assert result.returncode != 0
    assert "ENST00000334810" in result.stderr
    assert "not found" in result.stderr
