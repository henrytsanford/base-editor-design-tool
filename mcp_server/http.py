"""Serving the MCP server over HTTP, mounted on the web app.

The transport lives here rather than in server.py, so that module stays a list of what
the model is offered, and so the web app's whole dependency on the MCP stack is one
guarded import of one name.

Streamable HTTP, stateless, JSON responses. Stateless because Cloud Run runs this with
--min-instances 0 --max-instances 1, where a session id would bind a conversation to an
instance that can be gone before the next call. JSON rather than SSE because every one
of these tools answers with a single result, and a long-lived stream would have to
survive both Firebase Hosting in front and the host app's BaseHTTPMiddleware wrapping
it. What that costs is sampling, elicitation, roots and resumable streams -- none of
which the seven tools use.
"""
import logging

from . import tools
from .runs import RunStore
from .server import build_server

log = logging.getLogger(__name__)

# Where the sub-app is mounted, and where the endpoint sits inside it. '/api' plus
# '/mcp' makes '/api/mcp' the canonical URL, and Starlette 307s '/api/mcp/' onto it.
# Mounting at '/api/mcp' with the endpoint at '/' inverts that: the trailing slash
# becomes canonical and every client has to be told about it.
MOUNT_PATH = '/api'
MCP_PATH = '/mcp'

# 0.0.0.0, not the SDK's default '127.0.0.1'. That default auto-enables DNS-rebinding
# protection against an allow-list of localhost names, which answers 421 to every
# request carrying a real Host header -- including, measurably, a test client's. This
# is safe here for the reason `uvicorn --host 0.0.0.0` in the Dockerfile is safe: the
# container's network boundary is the boundary. Passing TransportSecuritySettings()
# instead would be worse than either: it defaults to protection on with empty
# allow-lists, and rejects everything.
MCP_HOST = '0.0.0.0'


def mount_mcp(app, settings):
    """Mounts the MCP server on `app`, returning the lifespan the host must enter.

    Starlette does not run a mounted sub-app's lifespan, and for this sub-app that
    lifespan is the task group its session manager owns -- without it the first POST
    answers 500, 'Task group is not initialized'. So the caller gets it back and is
    responsible for entering it.

    The context is built when that lifespan is entered, which is after app.state is
    populated, so the tools answer from the web app's own References, LocalStorage,
    ResultCache and JobPool rather than a second set of each.
    """
    def context():
        return tools.AppContext(
            settings=settings,
            references=app.state.references,
            runs=RunStore(settings, app.state.references,
                          storage=app.state.storage, tables=app.state.tables),
            pool=app.state.pool,
            limiter=app.state.limiter)

    mcp = build_server(settings, context=context)
    # Called exactly once per process: it builds the session manager the returned
    # lifespan runs, and that manager refuses a second run().
    sub = mcp.streamable_http_app(streamable_http_path=MCP_PATH, stateless_http=True,
                                  json_response=True, host=MCP_HOST)
    app.mount(MOUNT_PATH, sub)
    log.info('MCP over streamable HTTP at %s%s', MOUNT_PATH, MCP_PATH)
    return lambda: sub.router.lifespan_context(sub)
