"""Unit tests for Ensembl ID normalization and the retrying request helper.

These never touch the network: the module-level requests.Session is replaced
with a stub that replays a canned list of responses.
"""
import pytest
import requests

import base_editing_guide_designs as bed


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
    assert bed.strip_tr_version(raw) == expected


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
    monkeypatch.setattr(bed.time, "sleep", lambda seconds: None)


def install(monkeypatch, outcomes):
    session = FakeSession(outcomes)
    monkeypatch.setattr(bed, "_session", session)
    return session


def test_retries_then_succeeds(monkeypatch, no_sleep):
    session = install(monkeypatch, [FakeResponse(503), FakeResponse(200)])
    r = bed.ensembl_get("/lookup/id/ENST1", {})
    assert r.status_code == 200
    assert session.calls == 2


def test_gives_up_after_attempts(monkeypatch, no_sleep):
    session = install(monkeypatch, [FakeResponse(503)] * 5)
    with pytest.raises(bed.EnsemblUnavailable):
        bed.ensembl_get("/lookup/id/ENST1", {}, attempts=5)
    assert session.calls == 5


def test_4xx_returned_not_retried(monkeypatch, no_sleep):
    # A versioned ID returns 400; get_exons also relies on 4xx meaning
    # "no CDS mapping" rather than a transient failure.
    session = install(monkeypatch, [FakeResponse(400)])
    r = bed.ensembl_get("/lookup/id/ENST1.13", {})
    assert r.status_code == 400
    assert not r.ok
    assert session.calls == 1


def test_429_is_retried(monkeypatch, no_sleep):
    session = install(monkeypatch, [FakeResponse(429), FakeResponse(200)])
    assert bed.ensembl_get("/lookup/id/ENST1", {}).status_code == 200
    assert session.calls == 2


def test_retry_after_header_is_honoured(monkeypatch):
    slept = []
    monkeypatch.setattr(bed.time, "sleep", lambda seconds: slept.append(seconds))
    install(monkeypatch, [FakeResponse(429, {"Retry-After": "7"}), FakeResponse(200)])
    bed.ensembl_get("/lookup/id/ENST1", {})
    assert slept == [7.0]


def test_backoff_is_capped(monkeypatch):
    """An absurd Retry-After must not park the run for an hour."""
    slept = []
    monkeypatch.setattr(bed.time, "sleep", lambda seconds: slept.append(seconds))
    install(monkeypatch, [FakeResponse(503, {"Retry-After": "3600"}), FakeResponse(200)])
    bed.ensembl_get("/lookup/id/ENST1", {})
    assert slept == [bed.ENSEMBL_MAX_BACKOFF]


def test_exponential_backoff_is_capped(monkeypatch):
    slept = []
    monkeypatch.setattr(bed.time, "sleep", lambda seconds: slept.append(seconds))
    install(monkeypatch, [FakeResponse(503)] * 9)
    with pytest.raises(bed.EnsemblUnavailable):
        bed.ensembl_get("/lookup/id/ENST1", {}, attempts=9)
    assert all(s <= bed.ENSEMBL_MAX_BACKOFF for s in slept)
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
    assert bed.ensembl_get("/lookup/id/ENST1", {}).status_code == 200
    assert session.calls == 3


def test_a_timeout_is_always_passed(monkeypatch, no_sleep):
    """The pre-fix code passed no timeout, so a hung connection blocked forever."""
    seen = {}

    class RecordingSession:
        def get(self, url, headers=None, timeout=None):
            seen["timeout"] = timeout
            seen["url"] = url
            return FakeResponse(200)

    monkeypatch.setattr(bed, "_session", RecordingSession())
    bed.ensembl_get("/lookup/id/ENST1", {})
    assert seen["timeout"] == bed.ENSEMBL_TIMEOUT
    assert seen["url"] == "https://rest.ensembl.org/lookup/id/ENST1"
