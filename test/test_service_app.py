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
DESIGNS_URL = '/designs?transcript=%s&preset=ABE7.10' % ISY1


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


@pytest.fixture
def cached_url(client):
    """A finished design, for the tests that are about reading one back."""
    wait_for_table(client, DESIGNS_URL)
    return DESIGNS_URL


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


def test_the_second_request_is_served_from_the_cache(client, cached_url, tmp_path):
    """Requesting an uncached result twice runs it once (the M1 criterion)."""
    url = cached_url
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


# --- The search page, the table and downloads (slice 2) -------------------------


def test_the_search_page_submits_a_plain_get(client):
    """No script anywhere in the app, which is what keeps script-src 'none'."""
    response = client.get('/')
    assert response.status_code == 200
    assert 'action="/genes"' in response.text and 'method="get"' in response.text
    assert '<script' not in response.text


def test_the_preset_list_comes_from_the_engine(client):
    """The dropdown and the validator read one table, so they cannot drift."""
    from bedesign.engine import BE_TYPES
    body = client.get('/').text
    for name in BE_TYPES:
        assert 'value="%s"' % name in body


def test_the_form_preselects_the_default_editor(client):
    """An untouched form submits the same editor a bare /designs URL would run."""
    from bedesign.engine import DEFAULT_BE_TYPE
    body = client.get('/').text
    assert 'value="%s" selected' % DEFAULT_BE_TYPE in body
    assert body.count(' selected') == 1


def test_a_gene_prefix_lists_the_genes_it_could_mean(client):
    response = client.get('/genes?q=MAP2K')
    assert response.status_code == 200
    assert 'MAP2K1' in response.text and 'MAP2K7' in response.text


def test_a_gene_lists_its_transcripts_with_mane_first(client):
    response = client.get('/genes?q=MAP2K1&preset=ABE7.10')
    assert response.status_code == 200
    assert 'MANE Select' in response.text
    assert MAP2K1 in response.text
    # The link out carries the editor choice, so /designs needs nothing added.
    assert 'preset=ABE7.10' in response.text


def test_a_gene_search_that_matches_nothing_says_so(client):
    response = client.get('/genes?q=ZZZZZZZZ')
    assert response.status_code == 200
    assert 'No gene symbol starts with' in response.text


def test_a_hostile_gene_query_is_refused_and_escaped(client):
    response = client.get('/genes?q=%3Cscript%3E')
    assert response.status_code == 400
    assert '<script>' not in response.text


def test_no_page_in_the_app_needs_javascript(client):
    """script-src is 'none', so a page that needed a script would simply break."""
    urls = ['/', '/genes?q=MAP2K1', '/designs?transcript=%s&preset=ABE7.10' % ISY1]
    for url in urls:
        response = client.get(url)
        assert '<script' not in response.text, url
        assert "script-src 'none'" in response.headers['content-security-policy']


def test_filtering_a_cached_table_does_not_start_a_job(client, cached_url, tmp_path):
    """A filter is a view of a result, not a different result (design doc 3.1)."""
    url = cached_url
    manifests = list((tmp_path / 'results').rglob('manifest.json'))
    assert len(manifests) == 1

    started = client.app.state.pool.submit
    calls = []
    client.app.state.pool.submit = lambda *a, **k: calls.append(a) or started(*a, **k)
    filtered = client.get(url + '&mutation=Missense')
    assert filtered.status_code == 200
    assert 'guides match' in filtered.text
    assert calls == [], 'filtering must not start a job'
    assert len(list((tmp_path / 'results').rglob('manifest.json'))) == 1


def test_a_filter_the_result_cannot_satisfy_is_refused_not_silently_empty(
        client, cached_url):
    response = client.get(cached_url + '&mutation=Nonexistent')
    assert response.status_code == 400
    assert 'Nonexistent' in response.text


def test_a_header_link_sorts_the_table(client, cached_url):
    plain = _first_cell(client.get(cached_url).text)
    sorted_desc = _first_cell(
        client.get(cached_url + '&sort=%23+edits&dir=desc').text)
    assert plain and sorted_desc and plain != sorted_desc


def test_the_second_page_shows_different_guides(client, cached_url):
    assert 'Page 1 of' in client.get(cached_url).text
    assert (_first_cell(client.get(cached_url).text)
            != _first_cell(client.get(cached_url + '&page=2').text))


def test_a_download_returns_the_stored_file(client, cached_url):
    response = client.get('/designs/download?transcript=%s&preset=ABE7.10&file=designs'
                          % ISY1)
    assert response.status_code == 200
    assert response.headers['content-type'] == 'application/gzip'
    assert ISY1 in response.headers['content-disposition']
    assert 'attachment' in response.headers['content-disposition']
    assert response.content[:2] == b'\x1f\x8b'


def test_a_download_of_an_uncached_result_starts_no_job(client, tmp_path):
    """Otherwise a download would be a way around the limiter on /designs."""
    response = client.get('/designs/download?transcript=%s&preset=BE4max&file=designs'
                          % MAP2K1)
    assert response.status_code == 404
    assert not list((tmp_path / 'results').rglob('manifest.json'))


@pytest.mark.parametrize('query', [
    'file=../../etc/passwd',
    'file=manifest',
    'mutation=Missense',
])
def test_a_bad_download_request_is_refused(client, query):
    response = client.get('/designs/download?transcript=%s&preset=ABE7.10&%s'
                          % (ISY1, query))
    assert response.status_code == 400


def _first_cell(html):
    """The first data cell of the rendered table, for comparing orderings."""
    match = re.search(r'<tbody>\s*<tr><td>(.*?)</td>', html, re.S)
    return match.group(1) if match else ''
