"""Fixtures for the MCP server tests.

The server is exercised through a real in-process MCP client rather than by calling
the tool bodies: what is being checked is the surface a model actually sees, which
includes argument validation, the output schema and how an error comes back.

The client is driven synchronously. The SDK is async, but the rest of this suite is
not, and a blocking portal keeps these tests reading like their neighbours instead of
introducing a second testing style for one directory.
"""
import pytest

from service.config import Settings

try:
    import anyio.from_thread
    from mcp import Client

    from mcp_server.server import build_server
except ImportError:
    # The MCP stack is optional, the way the live Ensembl API and the reference
    # bundle are: it lives in requirements-mcp.txt, not requirements.txt. Skipping
    # collection rather than failing it keeps `pytest test/ -q` working for someone
    # who installed only the CLI's dependencies -- otherwise this conftest would
    # take the whole suite down, golden tests included.
    collect_ignore_glob = ['test_mcp_*.py']


@pytest.fixture(scope='session')
def portal():
    """One event loop, running on a thread, for every synchronous client call."""
    with anyio.from_thread.start_blocking_portal('asyncio') as running:
        yield running


@pytest.fixture(scope='module')
def mcp_settings(tmp_path_factory, refdata, clinvar_db):
    """Settings pointing at the test ClinVar database, not the bundle's.

    This matters more than it looks. The golden fixtures were produced against
    test/data/variant_summary_sample.txt.gz, so MAP2K1's 'Clinical significance'
    column only reproduces against that database. Pointed at the bundle's full
    ClinVar the designs still come out, with different annotations -- and KRTAP25-1
    happens to match either way, so a spot check would not notice.
    """
    return Settings(refdata=refdata, clinvar_db=str(clinvar_db),
                    results_dir=str(tmp_path_factory.mktemp('mcp_results')))


@pytest.fixture(scope='module')
def mcp_server(mcp_settings):
    return build_server(mcp_settings)


def text_of(result):
    """The text of a tool result's content parts.

    How an MCP result carries its text is an SDK detail, so it is unpacked in one
    place rather than in every test that reads an error message.
    """
    return ' '.join(getattr(part, 'text', '') for part in result.content)


class SyncClient(object):
    """A blocking view of the async MCP client.

    Only the calls these tests make. `raise_exceptions` is left off so a ToolError
    arrives as a result with `is_error` set -- which is what a model sees, and what
    the error-path tests are about.
    """

    def __init__(self, portal, client):
        self._portal = portal
        self._client = client

    def call_tool(self, name, arguments=None):
        return self._portal.call(self._client.call_tool, name, arguments or {})

    def call(self, name, arguments=None):
        """The structured payload of a call that is expected to succeed."""
        result = self.call_tool(name, arguments)
        assert not result.is_error, '%s failed: %s' % (name, text_of(result))
        return result.structured_content

    def error(self, name, arguments=None):
        """The message a call that is expected to fail shows the model."""
        result = self.call_tool(name, arguments)
        assert result.is_error, '%s unexpectedly succeeded' % name
        return text_of(result)

    def list_tools(self):
        return self._portal.call(self._client.list_tools)


@pytest.fixture
def client(mcp_server, portal):
    with portal.wrap_async_context_manager(Client(mcp_server)) as connected:
        yield SyncClient(portal, connected)
