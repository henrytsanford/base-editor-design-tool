"""The MCP server as the web app serves it: one URL, no session, designs in the pool.

These drive the mounted endpoint over HTTP rather than the in-process client that
test_mcp_tools.py uses, because what is under test is the front door -- routing, the
Host header, the rate limiter, the pool -- and none of that exists on stdio.
"""
import json
import os

import pytest
from starlette.testclient import TestClient

from service.app import create_app
from service.cachekey import MANIFEST
from service.config import Settings
from service.jobs import PoolBusy

# The app opens the reference bundle when it starts, so every test here needs one --
# the same reason test_service_app.py marks its whole module.
pytestmark = pytest.mark.bundle

MCP_URL = '/api/mcp'
HEADERS = {'Accept': 'application/json', 'Content-Type': 'application/json'}
ISY1 = 'ENST00000393295'
PRESET = 'ABE7.10'


@pytest.fixture
def settings(tmp_path, refdata, clinvar_db):
    return Settings(refdata=refdata, clinvar_db=str(clinvar_db),
                    results_dir=str(tmp_path / 'results'),
                    rate_burst=50, rate_seconds=1, job_timeout=120, mcp_wait=5)


@pytest.fixture
def client(settings):
    with TestClient(create_app(settings)) as connected:
        yield connected


def rpc(client, method, params=None, **kwargs):
    body = {'jsonrpc': '2.0', 'id': 1, 'method': method}
    if params is not None:
        body['params'] = params
    headers = dict(HEADERS, **kwargs.pop('headers', {}))
    return client.post(MCP_URL, json=body, headers=headers, **kwargs)


def call(client, name, arguments):
    """A tool call that is expected to succeed."""
    result = rpc(client, 'tools/call',
                 {'name': name, 'arguments': arguments}).json()['result']
    assert not result.get('isError'), result['content']
    return result['structuredContent']


def call_error(client, name, arguments):
    """The message a failing tool call shows the model."""
    result = rpc(client, 'tools/call',
                 {'name': name, 'arguments': arguments}).json()['result']
    assert result.get('isError'), '%s unexpectedly succeeded' % name
    return ' '.join(part.get('text', '') for part in result['content'])


# ---- the front door ----------------------------------------------------------

def test_the_endpoint_answers_at_one_url(client):
    listed = rpc(client, 'tools/list')
    assert listed.status_code == 200
    assert listed.headers['content-type'].startswith('application/json')
    assert len(listed.json()['result']['tools']) == 7
    # Starlette redirects the trailing slash rather than serving both, so the URL
    # clients are given is the one without it.
    redirected = client.post(MCP_URL + '/', json={'jsonrpc': '2.0', 'id': 1,
                                                  'method': 'tools/list'},
                             headers=HEADERS, follow_redirects=False)
    assert redirected.status_code == 307


def test_a_real_host_header_is_not_refused(client):
    """The regression test for the SDK's localhost default.

    streamable_http_app defaults to host='127.0.0.1', which auto-enables DNS-rebinding
    protection and answers 421 to every request carrying a real Host header -- which,
    behind Firebase Hosting, is every request there is.
    """
    answered = rpc(client, 'tools/list', headers={'Host': 'bedesigner.web.app'})
    assert answered.status_code == 200, answered.text


def test_no_session_is_handed_out(client):
    """Stateless, because Cloud Run scales this to zero at one instance: a session id
    would bind a conversation to an instance that can be gone before the next call."""
    assert 'mcp-session-id' not in rpc(client, 'tools/list').headers


def test_a_browser_preflight_is_not_answered(client):
    """CORS is deliberately absent: every intended client is server-side.

    Pinned so that adding CORS later is a deliberate act rather than a side effect.
    """
    assert client.options(MCP_URL).status_code == 405


def test_the_mcp_response_carries_the_security_headers(client):
    """The host app's middleware wraps the mount too."""
    assert rpc(client, 'tools/list').headers['x-frame-options'] == 'DENY'


def test_healthz_reports_the_mount(client):
    health = client.get('/healthz').json()
    assert health['mcp'] is True
    assert health['status'] == 'ok'


# ---- designs go through the pool ---------------------------------------------

def test_a_design_runs_in_the_job_pool(client, monkeypatch):
    """Not on an event-loop thread: the pool is what keeps CPU-bound design work out
    of the process serving the web app, and what enforces the wall-clock cap."""
    submitted = []
    real = client.app.state.pool.submit

    def record(key, target, params, client_=None):
        submitted.append(target)
        return real(key, target, params, client_)

    monkeypatch.setattr(client.app.state.pool, 'submit', record)
    run = call(client, 'design_guides', {'transcript_id': ISY1, 'preset': PRESET})
    assert [t.transcript_id for t in submitted] == [ISY1]
    assert run['counts']['designs'] > 0
    # The symbol the worker deliberately does not resolve, filled in off the index.
    assert run['target']['gene'] == 'ISY1'


def test_both_front_doors_are_one_result(client, settings):
    """The property the whole mount turns on.

    Both compute the same cache key, so a design asked for over MCP is the design the
    browser already has -- one job, one result on disk, no second wait.
    """
    call(client, 'design_guides', {'transcript_id': ISY1, 'preset': PRESET})
    page = client.get('/designs?transcript=%s&preset=%s' % (ISY1, PRESET))
    assert page.status_code == 200
    assert 'sgRNA sequence' in page.text, 'the browser was made to wait for a design'

    manifests = [os.path.join(root, MANIFEST)
                 for root, _dirs, files in os.walk(settings.results_dir)
                 if MANIFEST in files]
    assert len(manifests) == 1, 'two front doors wrote two results'


def test_the_worker_records_the_file_measurements(client, settings):
    """Measured where the uncompressed bytes already exist, so export_run does not
    decompress a whole file to report its size."""
    run = call(client, 'design_guides', {'transcript_id': ISY1, 'preset': PRESET})
    stored = [os.path.join(root, MANIFEST)
              for root, _dirs, files in os.walk(settings.results_dir)
              if MANIFEST in files]
    manifest = json.loads(open(stored[0]).read())
    assert manifest['stats']['designs']['sha256']
    exported = call(client, 'export_run',
                    {'run_id': run['run_id'], 'file': 'designs'})
    assert exported['sha256'] == manifest['stats']['designs']['sha256']
    assert exported['bytes'] == manifest['stats']['designs']['bytes']


def test_a_sequence_designs_through_the_pool(client):
    """The pool learned the nucleotide path, which only the MCP door can reach."""
    run = call(client, 'design_guides',
               {'sequence': 'ACGT' * 60, 'sequence_name': 'toy'})
    assert run['counts']['designs'] > 0
    assert run['files'] == ['designs', 'errors']
    assert run['clinvar'] == 'none'
    assert 'no gene' in call_error(client, 'export_run',
                                  {'run_id': run['run_id'], 'file': 'clinvar'})


# ---- what a busy or limited service says -------------------------------------

def test_the_rate_limit_turns_a_second_design_away(tmp_path, refdata, clinvar_db):
    """One token, so the second distinct design is refused -- and told what to do."""
    tight = Settings(refdata=refdata, clinvar_db=str(clinvar_db),
                     results_dir=str(tmp_path / 'r'), rate_burst=1,
                     rate_seconds=3600, job_timeout=120, mcp_wait=0)
    with TestClient(create_app(tight)) as client:
        call(client, 'design_guides', {'transcript_id': ISY1, 'preset': PRESET})
        message = call_error(client, 'design_guides',
                             {'transcript_id': ISY1, 'preset': 'ABE7.9'})
        assert 'Too many designs started from this address' in message
        assert 'not limited' in message, 'should say reading a run is still allowed'


def test_a_cached_run_is_never_charged(tmp_path, refdata, clinvar_db):
    """The same design twice costs one token, as a cache hit does on /designs."""
    tight = Settings(refdata=refdata, clinvar_db=str(clinvar_db),
                     results_dir=str(tmp_path / 'r'), rate_burst=1,
                     rate_seconds=3600, job_timeout=120, mcp_wait=0)
    with TestClient(create_app(tight)) as client:
        first = call(client, 'design_guides', {'transcript_id': ISY1, 'preset': PRESET})
        again = call(client, 'design_guides', {'transcript_id': ISY1, 'preset': PRESET})
        assert again['run_id'] == first['run_id']


def test_a_busy_pool_asks_the_model_to_call_again(tmp_path, refdata, clinvar_db,
                                                 monkeypatch):
    """No queue, by design. A browser retries through busy.html; a model has to be
    told to, in words it can act on.

    Its own app because mcp_wait has to be 0 -- Settings is frozen, and waiting the
    default 20 s to be told the pool is busy is not a thing to do in a test.
    """
    impatient = Settings(refdata=refdata, clinvar_db=str(clinvar_db),
                         results_dir=str(tmp_path / 'r'), rate_burst=50,
                         rate_seconds=1, mcp_wait=0)

    def always_busy(*args, **kwargs):
        raise PoolBusy()

    with TestClient(create_app(impatient)) as client:
        monkeypatch.setattr(client.app.state.pool, 'submit', always_busy)
        message = call_error(client, 'design_guides',
                             {'transcript_id': ISY1, 'preset': PRESET})
    assert 'Every design worker is busy' in message
    assert 'Nothing was lost' in message


def test_a_timed_out_design_names_what_to_narrow(client, monkeypatch):
    """The pool's own message ends 'Reload to try again', which a model cannot do."""
    monkeypatch.setattr(client.app.state.pool, 'failure', lambda key: 'stopped')
    monkeypatch.setattr(client.app.state.pool, 'failure_reason', lambda key: 'timeout')
    message = call_error(client, 'design_guides',
                         {'transcript_id': ISY1, 'preset': PRESET})
    assert 'one deaminase' in message
    assert 'Reload' not in message
