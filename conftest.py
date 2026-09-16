"""Shared test setup. Its location also puts the repository root on sys.path, so
test/ can import the script.

Markers: `ensembl` tests query the live REST API and run only with ENSEMBL_TESTS=1;
`bundle` tests need a local reference bundle under $REFDATA (default: refdata).
"""
import os
import re
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools'))

REFDATA_DIR = os.environ.get('REFDATA', 'refdata')

# Real ClinVar rows for the sample and golden-panel genes, in variant_summary.txt
# format: MAP2K1 659, HGS 152, PSMB5 42, KRTAP25-1 16, ISY1 12.
CLINVAR_FIXTURE = 'test/data/variant_summary_sample.txt.gz'


def pytest_addoption(parser):
    parser.addoption(
        '--regenerate-golden', action='store_true',
        help='rewrite the test/golden fixtures from the current code instead of '
             'comparing against them')


def pytest_configure(config):
    config.addinivalue_line('markers', 'ensembl: queries the live Ensembl REST API')
    config.addinivalue_line('markers', 'bundle: needs a local reference bundle')


def pytest_collection_modifyitems(config, items):
    import bedesign.transcript_source as transcript_source
    try:
        transcript_source.find_bundle(REFDATA_DIR)
        have_bundle = True
    except transcript_source.BundleNotFound:
        have_bundle = False
    for item in items:
        if item.get_closest_marker('ensembl') and not os.environ.get('ENSEMBL_TESTS'):
            item.add_marker(pytest.mark.skip(
                reason='queries the Ensembl REST API; set ENSEMBL_TESTS=1 to run'))
        if item.get_closest_marker('bundle') and not have_bundle:
            item.add_marker(pytest.mark.skip(
                reason='needs a reference bundle under %s; build one with '
                       'tools/build_reference.py' % REFDATA_DIR))


@pytest.fixture(scope='session')
def refdata():
    return REFDATA_DIR


@pytest.fixture(scope='session')
def clinvar_db(tmp_path_factory):
    """A ClinVar database over the sample genes, built by the real builder.

    Building it here rather than checking in a .db keeps the fixture honest: a
    change to build_clinvar.py that broke the schema or the column mapping would
    break these tests too.
    """
    import build_clinvar
    path = tmp_path_factory.mktemp('clinvar') / 'clinvar-test.db'
    build_clinvar.build(CLINVAR_FIXTURE, str(path))
    return str(path)


@pytest.fixture
def run_design(clinvar_db):
    """Runs the script, returning (result, output folder); removes the output afterwards.

    `--edit C-T` is a default that a caller can override by passing another
    `--edit` in `extra`; argparse keeps the last one.
    """
    folders = []

    def run(input_file, input_type, output_name, *extra, check=True):
        cmd = [sys.executable, 'base_editing_guide_designs.py',
               '--input-file', input_file, '--input-type', input_type,
               '--clinvar-db', clinvar_db, '--pam', 'NGG', '--intron-buffer', '30',
               '--edit', 'C-T', '--output-name', output_name, '--sg-len', '20', *extra]
        pattern = re.compile(re.escape(output_name) + r'_\d\d(-\d\d){5}$')
        # Only directories this call creates are ours to read or delete. A run
        # made by hand, or left behind by a failed test, can carry the same
        # name; picking it up would compare against, or delete, output this run
        # never wrote.
        before = {d for d in os.listdir('.') if pattern.match(d)}
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        found = sorted({d for d in os.listdir('.') if pattern.match(d)} - before)
        # registered before the assert, so a failing run is still cleaned up
        folders.extend(found)
        if check:
            assert result.returncode == 0, (
                f'Script failed with return code {result.returncode}\n'
                f'STDOUT: {result.stdout}\n'
                f'STDERR: {result.stderr}'
            )
        assert len(found) <= 1, f'run produced several output folders: {found}'
        return result, found[0] if found else None

    yield run
    for d in folders:
        shutil.rmtree(d, ignore_errors=True)
