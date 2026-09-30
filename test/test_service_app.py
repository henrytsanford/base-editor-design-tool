"""The /designs route, end to end against a real bundle.

Marked `bundle` so it skips where there is no reference data, like the golden tests.
ClinVar comes from the small `clinvar_db` fixture rather than the 767 MB database, and
ISY1 is the cheapest gene in the panel (0.6 s in the design doc's table).
"""
import html
import re
import time

import pytest

from bedesign import DESIGN_COLUMNS
from service.app import TABLE_COLUMNS, create_app
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


def test_health_answers_at_both_paths(client):
    # '/health' for anything outside, since Cloud Run never forwards '/healthz';
    # '/healthz' for a probe, which reaches the container directly and would restart
    # it on a 404. See the comment on the route in app.py.
    for path in ('/health', '/healthz'):
        assert client.get(path).json()['engine_version'], path


def test_health_reports_the_engine_version(client):
    body = client.get('/health').json()
    assert body['status'] == 'ok' and body['engine_version']


def test_robots_turns_crawlers_away(client):
    """A GET can start a job, so crawlers are refused here and by the meta tag."""
    body = client.get('/robots.txt').text
    assert 'Disallow: /' in body


def test_the_footer_names_the_reference_data(client):
    references = client.app.state.references
    body = client.get('/').text
    assert 'Ensembl release %s,' % references.release in body
    assert 'ClinVar %s.' % references.clinvar_version in body
    assert 'unknown' not in (references.release, references.clinvar_version)


def test_every_response_carries_the_security_headers(client):
    response = client.get('/health')
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
    assert re.search(r'<strong>[\d,]+ guides?</strong>', response.text)


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


def test_a_second_design_from_one_client_waits_without_spending_tokens(
        refdata, clinvar_db, tmp_path, monkeypatch):
    """One job per client. The waiting page retries, so it must not charge the
    client's bucket each time or the retries would rate-limit it."""
    from fastapi.testclient import TestClient
    settings = Settings(refdata=refdata, clinvar_db=clinvar_db,
                        results_dir=str(tmp_path / 'results'), max_jobs=2,
                        rate_burst=1, rate_seconds=3600)
    with TestClient(create_app(settings)) as client:
        monkeypatch.setattr(client.app.state.pool, 'client_busy', lambda c: True)
        for _ in range(3):
            response = client.get('/designs?transcript=%s&preset=ABE7.10' % MAP2K1)
            assert response.status_code == 429
            assert 'already have a design running' in response.text
            assert response.headers['retry-after']
        assert client.app.state.limiter.allow('testclient') is True


def test_health_fails_after_a_worker_crash(client, monkeypatch):
    monkeypatch.setattr(client.app.state.pool, 'healthy', lambda: False)
    response = client.get('/health')
    assert response.status_code == 503
    assert response.json()['status'] == 'worker crashed'


def test_a_result_past_the_table_limit_is_offered_as_downloads_only(
        refdata, clinvar_db, tmp_path):
    """Parsing a result with hundreds of thousands of guides would cost the web
    process hundreds of megabytes, so it is never parsed at all."""
    from fastapi.testclient import TestClient
    settings = Settings(refdata=refdata, clinvar_db=clinvar_db,
                        results_dir=str(tmp_path / 'results'), max_jobs=2,
                        rate_burst=50, rate_seconds=1, table_max_rows=10)
    with TestClient(create_app(settings)) as client:
        response = wait_for_table_or_downloads(client, DESIGNS_URL)
        assert response.status_code == 200
        assert 'too many to show as a table' in response.text
        assert 'sgRNA sequence' not in response.text
        assert client.app.state.tables.stats()['frames'] == 0
        link = re.search(r'href="(/designs/download\?[^"]+)"', response.text).group(1)
        download = client.get(link.replace('&amp;', '&'))
        assert download.status_code == 200 and download.content[:2] == b'\x1f\x8b'


def wait_for_table_or_downloads(client, url, seconds=90):
    for _ in range(seconds * 2):
        response = client.get(url)
        if 'Download the full result' in response.text or 'sgRNA sequence' in response.text:
            return response
        time.sleep(0.5)
    raise AssertionError('design did not finish within %ds' % seconds)


def test_a_result_evicted_mid_view_is_recomputed_not_a_500(client, cached_url,
                                                            monkeypatch):
    """The manifest was read, then the designs file was gone: eviction ran between."""
    storage = client.app.state.storage
    real_get = storage.get

    def evicted(key):
        if key.endswith('designs.tsv.gz'):
            raise KeyError(key)
        return real_get(key)

    monkeypatch.setattr(storage, 'get', evicted)
    client.app.state.tables = type(client.app.state.tables)()  # nothing parsed yet
    response = client.get(cached_url)
    assert response.status_code == 200
    assert 'http-equiv="refresh"' in response.text


def test_results_are_evicted_to_the_budget_after_each_job(refdata, clinvar_db, tmp_path):
    from fastapi.testclient import TestClient
    settings = Settings(refdata=refdata, clinvar_db=clinvar_db,
                        results_dir=str(tmp_path / 'results'), max_jobs=2,
                        rate_burst=50, rate_seconds=1, results_max_mb=1)
    with TestClient(create_app(settings)) as client:
        client.app.state.storage.evict = evicting = _Recorder(client.app.state.storage.evict)
        wait_for_table(client, DESIGNS_URL)
        assert evicting.calls == [('results', 'manifest.json', 2 ** 20)]


class _Recorder(object):
    def __init__(self, fn):
        self.fn, self.calls = fn, []

    def __call__(self, *args):
        self.calls.append(args)
        return self.fn(*args)


def test_the_downloads_are_described_above_the_table(client, cached_url):
    """Someone new to base editing should not have to scroll past the table, or
    guess from a filename, to find the file they want."""
    text = client.get(cached_url).text
    # A description only in a hover title is one nobody reads.
    visible = re.sub(r'title="[^"]*"', '', text)
    table = visible.index('<table')
    for label, filename, description in [
            ('all guides', 'designs.tsv.gz', 'Most people want this file'),
            ('Guides matching known variants', 'clinvar.tsv.gz', 'known human variant'),
            ('Design problems', 'errors.tsv.gz', 'Often empty')]:
        for piece in (label, filename, description):
            assert -1 < visible.find(piece) < table, piece
    assert '&amp;mdash;' not in text


def test_starting_jobs_too_fast_is_rate_limited(refdata, clinvar_db, tmp_path,
                                                monkeypatch):
    """The limiter guards the path that starts work, not cached reads."""
    from fastapi.testclient import TestClient
    settings = Settings(refdata=refdata, clinvar_db=clinvar_db,
                        results_dir=str(tmp_path / 'results'), max_jobs=2,
                        rate_burst=1, rate_seconds=3600)
    with TestClient(create_app(settings)) as client:
        # The first job is still running when the second request arrives, so the
        # one-job-per-client check would answer first; this test is about the bucket.
        monkeypatch.setattr(client.app.state.pool, 'client_busy', lambda c: False)
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


def test_the_table_leads_with_what_the_guide_does(client, cached_url):
    """The columns that hold one value for the whole result are in the heading, not
    repeated on every row; the download keeps them."""
    text = client.get(cached_url).text
    headers = [html.unescape(name) for name in
               re.findall(r'<th[^>]*><a [^>]*>(.*?)</a></th>', text)]
    assert headers == list(TABLE_COLUMNS)
    assert 'Ensembl Gene ID' not in headers and 'Genome assembly' not in headers
    row = re.search(r'<tbody>\s*<tr>(.*?)</tr>', text, re.S).group(1)
    assert row.count('<td>') == len(TABLE_COLUMNS)


def test_every_table_column_is_one_the_engine_writes():
    """A renamed engine column would otherwise drop out of the page unnoticed."""
    assert set(TABLE_COLUMNS) <= set(DESIGN_COLUMNS)


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


def _matched(html):
    """The (N, M) of the 'N of M guides match' the table prints when filtered."""
    match = re.search(r'([\d,]+) of ([\d,]+) guides match', html)
    return tuple(int(group.replace(',', '')) for group in match.groups()) if match else None


def test_the_table_page_carries_both_coverage_panels(client):
    response = wait_for_table(client, DESIGNS_URL)
    assert 'class="viz gene"' in response.text
    assert 'class="viz mtx"' in response.text
    # Between the downloads and the filter form, which is where the page puts it.
    assert (response.text.index('Download the full result')
            < response.text.index('Guide coverage')
            < response.text.index('class="filters"'))


def test_clicking_an_exon_filters_the_table_to_it(client):
    response = wait_for_table(client, DESIGNS_URL)
    link = re.search(r'href="(/designs\?[^"]*exon=2[^"]*)"', response.text)
    assert link, 'the map offers no exon link'
    filtered = client.get(html.unescape(link.group(1)))
    assert filtered.status_code == 200
    matched, total = _matched(filtered.text)
    assert 0 < matched < total
    # The selection is read, and undone, as a chip above the table.
    assert re.search(r'<a href="[^"]*" title="Remove this filter">Exon 2 ', filtered.text)


def test_clicking_a_square_filters_and_marks_the_map(client):
    """The number on the square is the number the table then shows."""
    response = wait_for_table(client, DESIGNS_URL)
    link = re.search(r'href="(/designs\?[^"]*sub=([A-Za-z]{3}-[A-Za-z]{3})[^"]*)"',
                     response.text)
    assert link, 'the matrix offers no substitution link'
    url = html.unescape(link.group(1))
    promised = int(re.search(
        r'>[A-Za-z]+ to [A-Za-z]+ · (\d+) guide', response.text).group(1))
    filtered = client.get(url)
    assert filtered.status_code == 200
    assert _matched(filtered.text)[0] == promised
    # and the map above it fades every bin without one of those guides
    assert 'class="bin dim"' in filtered.text
    assert 'class="bin"' in filtered.text


def test_a_bar_on_the_map_links_to_its_exon(client):
    response = wait_for_table(client, DESIGNS_URL)
    bins = re.findall(r'<a href="([^"]*)"><g class="bin', response.text)
    assert bins, 'the map offers no bin link'
    assert all('exon=' in html.unescape(url) for url in bins)


def test_an_exon_and_a_square_combine(client):
    """The two chart selections narrow each other instead of replacing each other,
    and every link the charts offer then still matches a guide."""
    response = wait_for_table(client, DESIGNS_URL)
    sub = re.search(r'sub=([A-Za-z]{3}-[A-Za-z]{3})', response.text).group(1)
    selected = client.get(DESIGNS_URL + '&sub=' + sub)
    exon = re.search(r'<a href="([^"]*)"><g class="exon"', selected.text)
    assert exon, 'no exon holds the selected change'
    both = client.get(html.unescape(exon.group(1)))
    assert both.status_code == 200
    assert 'title="Remove this filter">Exon ' in both.text
    assert 'sub=' + sub in html.unescape(exon.group(1))
    assert _matched(both.text)[0] > 0


def test_a_square_matching_nothing_in_the_selected_exon_is_not_a_link(client):
    wait_for_table(client, DESIGNS_URL)
    response = client.get(DESIGNS_URL + '&exon=2')
    cells = re.findall(r'(<a href="[^"]*">)?<g class="cell( empty)?">', response.text)
    assert any(empty for _, empty in cells)
    assert not any(link and empty for link, empty in cells)


def test_links_land_where_the_click_is_answered(client):
    """A reload lands on the map for a chart click and on the filters for anything
    done to the table, rather than at the top of the page."""
    response = wait_for_table(client, DESIGNS_URL)
    text = response.text
    assert 'id="coverage"' in text and 'id="results"' in text
    assert re.search(r'href="/designs\?[^"]*exon=\d+#coverage"', text)
    assert re.search(r'href="/designs\?[^"]*sub=[A-Za-z-]+#coverage"', text)
    assert re.search(r'href="/designs\?[^"]*sort=[^"]*#results"', text)
    assert 'action="/designs#results"' in text
    # The downloads are files, not views, and carry no fragment.
    assert not re.search(r'href="/designs/download[^"]*#', text)


def test_a_substitution_this_editor_cannot_make_is_a_400(client):
    """Tryptophan to tryptophan needs no edit at all, so no guide makes it."""
    wait_for_table(client, DESIGNS_URL)
    response = client.get(DESIGNS_URL + '&sub=Trp-Trp')
    assert response.status_code == 400


@pytest.mark.parametrize('query', ['exon=0', 'exon=99999', 'sub=Xyz-Gly',
                                   'sub=drop+table', 'exon=two'])
def test_a_bad_coverage_filter_is_refused(client, query):
    response = client.get('%s&%s' % (DESIGNS_URL, query))
    assert response.status_code == 400


def test_the_filter_form_keeps_the_chart_selection(client):
    """Applying a filter narrows what was clicked instead of discarding it."""
    wait_for_table(client, DESIGNS_URL)
    response = client.get(DESIGNS_URL + '&exon=2')
    form = response.text[response.text.index('class="filters"'):
                         response.text.index('</form>')]
    assert '<input type="hidden" name="exon" value="2">' in form
    assert 'name="sub"' not in form


def test_a_chip_removes_only_its_own_filter(client):
    wait_for_table(client, DESIGNS_URL)
    response = client.get(DESIGNS_URL + '&exon=2&mutation=Missense')
    chips = dict((label.strip(), html.unescape(url)) for url, label in re.findall(
        r'<a href="([^"]*)" title="Remove this filter">([^<]*)<', response.text))
    category = 'Mutation category: Missense'
    assert set(chips) == {'Exon 2', category}
    assert 'exon=' not in chips['Exon 2'] and 'mutation=Missense' in chips['Exon 2']
    assert 'exon=2' in chips[category] and 'mutation=' not in chips[category]
    assert 'Clear all' in response.text


def test_the_map_key_is_the_consequence_filter(client):
    """Each class in the key links to exactly as many guides as it says it counts,
    and there is no second list of consequences in the form. Guides with no edit
    are not drawn, so theirs is the one entry without a swatch."""
    response = wait_for_table(client, DESIGNS_URL)
    assert 'name="mutation"' not in response.text
    assert 'c-none' not in response.text
    key = re.findall(r'<a href="([^"]*consequence=(\w+)[^"]*)"[^>]*>'
                     r'(?:<i class="c-\w+"></i>)?([\d,]+) ', response.text)
    assert {name for _, name, _ in key} >= {'mis', 'none'}
    for url, name, count in key:
        filtered = client.get(html.unescape(url))
        assert re.search(r'<strong>%s of [\d,]+ guides match' % count, filtered.text), name
        assert 'class="on" aria-current="true"' in filtered.text
        assert '<input type="hidden" name="consequence" value="%s">' % name \
            in filtered.text


def test_a_mutation_link_still_filters_and_survives_the_form(client):
    """The form no longer offers the engine's categories, but a link naming one
    keeps working, and applying another filter does not drop it."""
    wait_for_table(client, DESIGNS_URL)
    response = client.get(DESIGNS_URL + '&mutation=Missense')
    assert 'guides match' in response.text
    assert '<input type="hidden" name="mutation" value="Missense">' in response.text


def test_a_single_edit_editor_offers_no_edit_dropdown(client):
    """ABE7.10 makes only A-G, and a dropdown with one choice is not a choice."""
    response = wait_for_table(client, DESIGNS_URL)
    assert 'name="deaminase"' not in response.text


def test_selecting_on_a_chart_leaves_its_caption_alone(client):
    """A selection that rewrote the caption would move the drawing under it."""
    def captions(text):
        # The consequence links carry the selection in their URLs, which the reader
        # never sees, so only what is drawn is compared.
        return [re.sub(r' href="[^"]*"', '', caption)
                for caption in re.findall(r'<figcaption>.*?</figcaption>', text, re.S)]

    response = wait_for_table(client, DESIGNS_URL)
    sub = re.search(r'sub=([A-Za-z]{3}-[A-Za-z]{3})', response.text).group(1)
    for query in ('&exon=2', '&sub=' + sub):
        selected = client.get(DESIGNS_URL + query)
        # The map's key differs only in the one entry kept hidden but laid out.
        assert (captions(selected.text.replace(' class="off" aria-hidden="true"', ''))
                == captions(response.text.replace(' class="off" aria-hidden="true"', '')))
    assert re.search(r'change is made at \d+\s+site', selected.text)


def test_a_selected_square_shows_as_a_chip_in_three_letter_names(client):
    response = wait_for_table(client, DESIGNS_URL)
    sub = re.search(r'sub=([A-Za-z]{3}-[A-Za-z]{3})', response.text).group(1)
    selected = client.get(DESIGNS_URL + '&sub=' + sub)
    source, target = sub.split('-')
    label = '%s → %s' % (source, 'Stop' if target == 'Ter' else target)
    assert 'title="Remove this filter">%s ' % label in selected.text


def test_the_matrix_is_folded_until_a_square_is_selected(client):
    """Open whenever a square is selected, since that is what the map then shows."""
    response = wait_for_table(client, DESIGNS_URL)
    assert '<details class="panel matrix">' in response.text
    sub = re.search(r'sub=([A-Za-z]{3}-[A-Za-z]{3})', response.text).group(1)
    selected = client.get(DESIGNS_URL + '&sub=' + sub)
    assert '<details class="panel matrix" open>' in selected.text
