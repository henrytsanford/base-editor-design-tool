"""The MCP mount, and that the web app does not depend on it.

Deliberately not named test_mcp_*, so test/conftest.py does not skip this file when
the MCP stack is absent -- the whole point of the first test is to run in exactly that
case.
"""
import os
import subprocess
import sys

# The repo root, which is not the same path inside the container image as it is here.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MOUNT_PATH = '/api'

# Makes `import mcp` invisible, the way a requirements.txt-only install leaves it.
# In a subprocess because the guard in service/app.py runs when the module is first
# imported, and this suite has imported it already.
WITHOUT_MCP = '''
import importlib.util
real = importlib.util.find_spec
importlib.util.find_spec = (
    lambda name, *a, **k: None if name == 'mcp' else real(name, *a, **k))
import service.app
assert service.app.mount_mcp is None, 'the guard did not see mcp as missing'
app = service.app.create_app()
mounted = [r for r in app.routes if getattr(r, 'path', '') == %r]
assert not mounted, 'the app mounted an MCP route without the MCP stack'
print('ok')
''' % MOUNT_PATH


def test_the_app_imports_and_serves_without_the_mcp_stack():
    """requirements-mcp.txt is optional, so service/app.py must not need it.

    An unguarded import here would break `pytest test/ -q` and the CLI install for
    everyone, not just the MCP tests -- service/app.py is imported by the whole
    service suite.
    """
    done = subprocess.run([sys.executable, '-c', WITHOUT_MCP],
                          capture_output=True, text=True, cwd=ROOT)
    assert done.returncode == 0, done.stderr
    assert 'ok' in done.stdout


def test_the_mount_is_there_when_the_stack_is():
    """The other half: with mcp importable, the route exists.

    create_app opens no reference data -- References is built in the lifespan -- so
    this needs no bundle.
    """
    import service.app
    if service.app.mount_mcp is None:
        import pytest
        pytest.skip('mcp is not installed')
    app = service.app.create_app()
    assert any(getattr(r, 'path', '') == MOUNT_PATH for r in app.routes)
