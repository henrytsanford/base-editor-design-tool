"""Shared test setup. Its location also puts the repository root on sys.path, so
test/ can import the script.

Markers: `ensembl` tests query the live REST API and run only with ENSEMBL_TESTS=1;
`bundle` tests need a local reference bundle under $REFDATA (default: refdata).
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools'))

REFDATA_DIR = os.environ.get('REFDATA', 'refdata')


def pytest_configure(config):
    config.addinivalue_line('markers', 'ensembl: queries the live Ensembl REST API')
    config.addinivalue_line('markers', 'bundle: needs a local reference bundle')


def pytest_collection_modifyitems(config, items):
    import transcript_source
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
