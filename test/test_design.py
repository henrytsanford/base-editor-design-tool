import subprocess
import shutil
import os
import pandas as pd
from base_editing_guide_designs import get_aa_map

SAMPLE_DATA_DIR = "Sample_data"


def format_args(input_file, input_type, output_name):
    cmd = [
        "python",
        "base_editing_guide_designs.py",
        "--input-file",
        input_file,
        "--input-type",
        input_type,
        "--variant-file",
        "variant_summary.txt",
        "--pam",
        "NGG",
        "--intron-buffer",
        "30",
        "--edit",
        "C-T",
        "--output-name",
        output_name,
        "--sg-len",
        "20",
    ]

    return cmd


def assert_output_files_equal(output_dir, sample_data_results_dir, output_name):
    output_file = f"sgrna_designs_{output_name}"
    with open(f"{output_dir}/{output_file}.txt") as f1, open(
        f"{sample_data_results_dir}/{output_file}.txt"
    ) as f2:
        # don't test columns with known differences:
        # Ensembl Gene ID missing from test file
        # Clinical significance broken by ClinVar changes
        # BsmBI flag logic changed from test file
        f1_contents = (
            pd.read_csv(f1, sep="\t")
            .drop(columns=["Ensembl Gene ID", "Clinical significance", "BsmBI flag"])
            .to_csv(sep="\t", index=False)
        )
        f2_contents = (
            pd.read_csv(f2, sep="\t")
            .drop(columns=["Clinical significance", "BsmBI flag"])
            .to_csv(sep="\t", index=False)
        )

        # test file uses single character amino acid codes
        aa_map = get_aa_map()
        aa_map["Ter"] = "*"
        for key, value in aa_map.items():
            f1_contents = f1_contents.replace(key, value)

        assert f1_contents == f2_contents
    with open(f"{output_dir}/README.txt") as f1, open(
        f"{sample_data_results_dir}/README.txt"
    ) as f2:
        # final two lines contain timestamps, don't check them
        assert f1.readlines()[:5] == f2.readlines()[:5]


def test_fasta_input():
    output_name = "GFP"
    cmd = format_args(f"{SAMPLE_DATA_DIR}/GFP.fasta", "nuc", output_name)
    output_dir = ""
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60,
        )

        assert result.returncode == 0, (
            f"Script failed with return code {result.returncode}\n"
            f"STDOUT: {result.stdout}\n"
            f"STDERR: {result.stderr}"
        )

        output_dir = [c for c in os.listdir(".") if c.startswith(output_name)][0]
        sample_data_results_dir = f"{SAMPLE_DATA_DIR}/GFP_19-05-07-14-41-15"
        assert_output_files_equal(output_dir, sample_data_results_dir, output_name)

    finally:
        # clean up files
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)


def test_tid_input():
    output_name = "sample"
    cmd = format_args(f"{SAMPLE_DATA_DIR}/sample_input.txt", "tid", output_name)
    output_dir = ""
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )

        assert result.returncode == 0, (
            f"Script failed with return code {result.returncode}\n"
            f"STDOUT: {result.stdout}\n"
            f"STDERR: {result.stderr}"
        )

        output_dir = [c for c in os.listdir(".") if c.startswith(output_name)][0]
        sample_data_results_dir = f"{SAMPLE_DATA_DIR}/sample_19-05-07-14-35-42"
        assert_output_files_equal(output_dir, sample_data_results_dir, output_name)

    finally:
        # clean up files
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)