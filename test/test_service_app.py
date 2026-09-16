"""The /designs route, end to end against a real bundle.

Marked `bundle` so it skips where there is no reference data, like the golden tests.
ClinVar comes from the small `clinvar_db` fixture rather than the 767 MB database, and
ISY1 is the cheapest gene in the panel (0.6 s in the design doc's table).
"""
import re
import time

import pytest

from service.app import create_app
from service.config import Settings
from service.jobs import PoolBusy

pytestmark = pytest.mark.bundle

ISY1 = 'ENST00000393295'
MAP2K1 = 'ENST00000307102'


@pytest.fixture
def client(refdata, clinvar_db, tmp_path):
    from fastapi.testclient import TestClient
    settings = Settings(refdata=refdata, clinvar_db=clinvar_db,
                        results_dir=str(tmp_path / 'results'), max_jobs=2,
                        rate_burst=50, rate_seconds=1)
    with TestClient(create_app(settings)) as client:
        yield client


def wait_for_table(client, url, seconds=90):
    """Follows the running page the way a browser's meta refresh would."""
    for _ in range(seconds * 2):
        response = client.get(url)
        if 'sgRNA sequence' in response.text or response.status_code >= 400:
            return response
        time.sleep(0.5)
    raise AssertionError('design did not finish within %ds' % seconds)


def test_healthz_reports_the_engine_version(client):
    body = client.get('/healthz').json()
    assert body['status'] == 'ok' and body['engine_version']


def test_robots_turns_crawlers_away(client):
    """A GET can start a job, so crawlers are refused here and by the meta tag."""
    body = client.get('/robots.txt').text
    assert 'Disallow: /' in body


def test_every_response_carries_the_security_headers(client):
    response = client.get('/healthz')
    assert response.headers['x-content-type-options'] == 'nosniff'
    assert response.headers['x-frame-options'] == 'DENY'
    assert response.headers['referrer-policy'] == 'no-referrer'
    policy = response.headers['content-security-policy']
    assert "default-src 'self'" in policy and "frame-ancestors 'none'" in policy
    assert 'unsafe-inline' not in policy and 'unsafe-eval' not in policy


def test_the_running_page_refreshes_without_javascript(client):
    """script-src is 'none', so the page must not need a script to re-check."""
    response = client.get('/designs?transcript=%s&preset=ABE7.10' % ISY1)
    assert response.status_code == 200
    assert 'http-equiv="refresh"' in response.text
    assert '<script' not in response.text


def test_a_design_runs_and_renders_its_guides(client):
    url = '/designs?transcript=%s&preset=ABE7.10' % ISY1
    response = wait_for_table(client, url)
    assert response.status_code == 200
    assert 'sgRNA sequence' in response.text
    assert re.search(r'\d+ guides?\.', response.text)


def test_the_second_request_is_served_from_the_cache(client, tmp_path):
    """Requesting an uncached result twice runs it once (the M1 criterion)."""
    url = '/designs?transcript=%s&preset=ABE7.10' % ISY1
    wait_for_table(client, url)
    manifests = list((tmp_path / 'results').rglob('manifest.json'))
    assert len(manifests) == 1
    first = manifests[0].read_text()

    started = client.app.state.pool.submit
    calls = []
    client.app.state.pool.submit = lambda *a, **k: calls.append(a) or started(*a, **k)
    assert 'sgRNA sequence' in client.get(url).text
    assert calls == [], 'a cached result should not start a job'
    assert manifests[0].read_text() == first


def test_a_different_parameter_is_a_different_result(client, tmp_path):
    wait_for_table(client, '/designs?transcript=%s&preset=ABE7.10' % ISY1)
    wait_for_table(client, '/designs?transcript=%s&preset=BE4max' % ISY1)
    assert len(list((tmp_path / 'results').rglob('manifest.json'))) == 2


def test_a_preset_and_its_explicit_parameters_share_one_result(client, tmp_path):
    """Design doc 3.1: the key is the resolved parameters, not how they were spelled."""
    wait_for_table(client, '/designs?transcript=%s&preset=ABE7.10' % ISY1)
    wait_for_table(client, '/designs?transcript=%s&pam=NGG&window=4-7&sg_len=20&edit=A-G'
                   % ISY1)
    assert len(list((tmp_path / 'results').rglob('manifest.json'))) == 1


@pytest.mark.parametrize('query', [
    'transcript=../../etc/passwd',
    'transcript=%s&pam=%%27+OR+1%%3D1--' % ISY1,
    'transcript=%s&preset=NOPE' % ISY1,
    'transcript=%s&sg_len=9999' % ISY1,
    'transcript=%s&evil=1' % ISY1,
    'transcript=ENST99999999999',
    '',
])
def test_a_bad_request_is_refused_without_starting_a_job(client, tmp_path, query):
    response = client.get('/designs?' + query)
    assert response.status_code == 400
    assert not list((tmp_path / 'results').rglob('manifest.json'))


def test_an_error_page_shows_no_traceback_or_path(client):
    response = client.get('/designs?transcript=ENST99999999999')
    assert response.status_code == 400
    assert 'Traceback' not in response.text
    assert 'refdata' not in response.text and '/Users' not in response.text


def test_an_unknown_transcript_is_refused_by_the_bundle(client):
    response = client.get('/designs?transcript=ENST99999999999')
    assert 'not in the reference bundle' in response.text


def test_html_from_the_data_is_escaped(client):
    """Autoescape is on; a value that looks like markup must arrive as text."""
    response = client.get('/designs?transcript=%s&preset=%s'
                          % (ISY1, '<script>alert(1)</script>'))
    assert response.status_code == 400
    assert '<script>' not in response.text
    assert '&lt;script&gt;' in response.text


def test_a_busy_pool_says_so_and_asks_for_a_retry(client, monkeypatch):
    def busy(*args, **kwargs):
        raise PoolBusy()

    monkeypatch.setattr(client.app.state.pool, 'submit', busy)
    response = client.get('/designs?transcript=%s&preset=ABE7.10' % MAP2K1)
    assert response.status_code == 503
    assert response.headers['retry-after']
    assert 'Busy' in response.text


def test_starting_jobs_too_fast_is_rate_limited(refdata, clinvar_db, tmp_path):
    """The limiter guards the path that starts work, not cached reads."""
    from fastapi.testclient import TestClient
    settings = Settings(refdata=refdata, clinvar_db=clinvar_db,
                        results_dir=str(tmp_path / 'results'), max_jobs=2,
                        rate_burst=1, rate_seconds=3600)
    with TestClient(create_app(settings)) as client:
        assert client.get('/designs?transcript=%s&preset=ABE7.10' % ISY1).status_code == 200
        second = client.get('/designs?transcript=%s&preset=BE4max' % ISY1)
        assert second.status_code == 429
        assert 'Too many designs' in second.text


def test_a_cached_result_is_not_rate_limited(refdata, clinvar_db, tmp_path):
    from fastapi.testclient import TestClient
    settings = Settings(refdata=refdata, clinvar_db=clinvar_db,
                        results_dir=str(tmp_path / 'results'), max_jobs=2,
                        rate_burst=1, rate_seconds=3600)
    with TestClient(create_app(settings)) as client:
        url = '/designs?transcript=%s&preset=ABE7.10' % ISY1
        wait_for_table(client, url)
        for _ in range(5):
            assert client.get(url).status_code == 200


def test_the_manifest_records_what_made_the_result(client, tmp_path):
    import json
    wait_for_table(client, '/designs?transcript=%s&preset=ABE7.10' % ISY1)
    manifest = json.loads(
        list((tmp_path / 'results').rglob('manifest.json'))[0].read_text())
    assert manifest['transcript_id'] == ISY1
    assert manifest['engine_version'] and manifest['reference']
    assert manifest['params']['edit'] == 'A-G'
    assert manifest['designs'] >= 0 and 'runtime_seconds' in manifest


def test_results_are_written_under_the_digest_only(client, tmp_path):
    """No user-supplied string is ever a path segment."""
    wait_for_table(client, '/designs?transcript=%s&preset=ABE7.10' % ISY1)
    manifest = list((tmp_path / 'results').rglob('manifest.json'))[0]
    parts = manifest.relative_to(tmp_path / 'results').parts
    assert parts[0] == 'results' and len(parts[2]) == 64
    assert all(re.fullmatch(r'[0-9a-f]+', p) for p in parts[1:3])
