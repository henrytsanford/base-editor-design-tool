"""The bundle lookups the app does before it schedules anything."""
import threading

import pytest

from service.references import References

pytestmark = pytest.mark.bundle

MAP2K1 = 'ENST00000307102'


@pytest.fixture
def references(refdata, clinvar_db):
    return References(refdata, clinvar_db)


def test_it_reports_the_release_and_clinvar_version(references):
    assert references.release
    assert references.clinvar_version


def test_a_known_transcript_exists(references):
    assert references.transcript_exists(MAP2K1) is True


def test_an_unknown_transcript_does_not(references):
    assert references.transcript_exists('ENST99999999999') is False


def test_lookups_are_correct_under_concurrency(references):
    """One sqlite3 connection shared across the ASGI thread pool interleaves its
    cursors, and a transcript that is in the bundle reads back as missing. Each
    thread gets its own connection, so this must be unanimous."""
    results = []
    lock = threading.Lock()

    def check():
        found = references.transcript_exists(MAP2K1)
        with lock:
            results.append(found)

    threads = [threading.Thread(target=check) for _ in range(100)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert results == [True] * 100
