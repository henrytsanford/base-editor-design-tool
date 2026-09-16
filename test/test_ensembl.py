"""Unit tests for Ensembl ID normalization and the retrying request helper.

These never touch the network: the module-level requests.Session in
transcript_source is replaced with a stub that replays a canned list of
responses.
"""
import pytest
import requests

from bedesign import engine
import bedesign.transcript_source as ts


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("ENST00000294952.13", "ENST00000294952"),
        ("ENST00000294952", "ENST00000294952"),
        ("ENSG00000164062.15", "ENSG00000164062"),
        ("ENSP00000294952.8", "ENSP00000294952"),
        ("ENSMUST00000000001.5", "ENSMUST00000000001"),
        # spreadsheet exports leave whitespace behind
        (" ENST00000294952.13\r", "ENST00000294952"),
        # not an Ensembl ID: left alone rather than mangled
        ("NM_001291281.2", "NM_001291281.2"),
    ],
)
def test_strip_tr_version(raw, expected):
    assert engine.strip_tr_version(raw) == expected


class FakeResponse:
    def __init__(self, status_code, headers=None):
        self.status_code = status_code
        self.headers = headers or {}
        self.text = "ACGT"

    @property
    def ok(self):
        return self.status_code < 400


class FakeSession:
    """Replays `outcomes` in order; each is a FakeResponse or an exception."""

    def __init__(self, outcomes):
        self.outcomes = list(outcomes)
        self.calls = 0

    def get(self, url, headers=None, timeout=None):
        self.calls += 1
        outcome = self.outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome


@pytest.fixture
def no_sleep(monkeypatch):
    monkeypatch.setattr(ts.time, "sleep", lambda seconds: None)


def install(monkeypatch, outcomes):
    session = FakeSession(outcomes)
    monkeypatch.setattr(ts, "_session", session)
    return session


def test_retries_then_succeeds(monkeypatch, no_sleep):
    session = install(monkeypatch, [FakeResponse(503), FakeResponse(200)])
    r = ts.ensembl_get("/lookup/id/ENST1", {})
    assert r.status_code == 200
    assert session.calls == 2


def test_gives_up_after_attempts(monkeypatch, no_sleep):
    session = install(monkeypatch, [FakeResponse(503)] * 5)
    with pytest.raises(ts.EnsemblUnavailable):
        ts.ensembl_get("/lookup/id/ENST1", {}, attempts=5)
    assert session.calls == 5


def test_4xx_returned_not_retried(monkeypatch, no_sleep):
    # A versioned ID returns 400; get_exons also relies on 4xx meaning
    # "no CDS mapping" rather than a transient failure.
    session = install(monkeypatch, [FakeResponse(400)])
    r = ts.ensembl_get("/lookup/id/ENST1.13", {})
    assert r.status_code == 400
    assert not r.ok
    assert session.calls == 1


def test_429_is_retried(monkeypatch, no_sleep):
    session = install(monkeypatch, [FakeResponse(429), FakeResponse(200)])
    assert ts.ensembl_get("/lookup/id/ENST1", {}).status_code == 200
    assert session.calls == 2


@pytest.mark.parametrize("retry_after,expected", [("7", 7.0), ("3600", ts.ENSEMBL_MAX_BACKOFF)])
def test_retry_after_is_honoured_up_to_the_cap(monkeypatch, retry_after, expected):
    """An absurd Retry-After must not park the run for an hour."""
    slept = []
    monkeypatch.setattr(ts.time, "sleep", slept.append)
    install(monkeypatch, [FakeResponse(429, {"Retry-After": retry_after}), FakeResponse(200)])
    ts.ensembl_get("/lookup/id/ENST1", {})
    assert slept == [expected]


def test_exponential_backoff_is_capped(monkeypatch):
    slept = []
    monkeypatch.setattr(ts.time, "sleep", lambda seconds: slept.append(seconds))
    install(monkeypatch, [FakeResponse(503)] * 9)
    with pytest.raises(ts.EnsemblUnavailable):
        ts.ensembl_get("/lookup/id/ENST1", {}, attempts=9)
    assert all(s <= ts.ENSEMBL_MAX_BACKOFF for s in slept)
    assert slept == sorted(slept)  # monotonically backing off


def test_timeouts_and_connection_errors_are_retried(monkeypatch, no_sleep):
    session = install(
        monkeypatch,
        [
            requests.exceptions.ConnectionError("reset"),
            requests.exceptions.Timeout("slow"),
            FakeResponse(200),
        ],
    )
    assert ts.ensembl_get("/lookup/id/ENST1", {}).status_code == 200
    assert session.calls == 3


def test_a_timeout_is_always_passed(monkeypatch, no_sleep):
    """The pre-fix code passed no timeout, so a hung connection blocked forever."""
    seen = {}

    class RecordingSession:
        def get(self, url, headers=None, timeout=None):
            seen["timeout"] = timeout
            seen["url"] = url
            return FakeResponse(200)

    monkeypatch.setattr(ts, "_session", RecordingSession())
    ts.ensembl_get("/lookup/id/ENST1", {})
    assert seen["timeout"] == ts.ENSEMBL_TIMEOUT
    assert seen["url"] == "https://rest.ensembl.org/lookup/id/ENST1"


class TestNonCodingTranscripts:
    """/map/cds returns 500 for a transcript with no CDS.

    Verified against the live API: ENST00000508832 (MALAT1) and ENST00000610481
    both return a 500 HTML page while a coding transcript returns 200 in the
    same window. Retrying that to exhaustion and raising would abort a whole run
    over one non-coding transcript, so it has to be told apart from an outage.
    """

    def test_persistent_500_means_no_cds_when_the_api_answers(self, monkeypatch, no_sleep):
        # seven 500s for /map/cds, then a healthy /lookup proving the API is up
        session = install(monkeypatch, [FakeResponse(500)] * 7 + [FakeResponse(200)])
        assert ts.EnsemblRestSource().cds_mappings("ENST1", 1000) is None
        assert session.calls == 8

    def test_persistent_500_still_raises_when_the_api_is_down(self, monkeypatch, no_sleep):
        # /map/cds fails, and so does the follow-up /lookup: a real outage
        install(monkeypatch, [FakeResponse(500)] * 7 + [FakeResponse(503)])
        with pytest.raises(ts.EnsemblUnavailable):
            ts.EnsemblRestSource().cds_mappings("ENST1", 1000)

    def test_4xx_still_means_no_cds(self, monkeypatch, no_sleep):
        session = install(monkeypatch, [FakeResponse(400)])
        assert ts.EnsemblRestSource().cds_mappings("ENST1", 1000) is None
        assert session.calls == 1

    def test_reachable_does_not_retry(self, monkeypatch, no_sleep):
        session = install(monkeypatch, [FakeResponse(503)])
        assert ts.EnsemblRestSource().reachable("ENST1") is False
        assert session.calls == 1
