"""The HTTP layer.

Design doc 3: the URL is the request. `GET /designs` validates, computes the cache
key, and then either renders the cached table, joins a running job, or starts one.
Nothing is stored about a request beyond the in-memory job table, so a restart costs
a recomputation rather than a broken link.

Slice 1 serves the design doc's 3 end to end. The search page, gene lookup, table
filtering and downloads are slice 2.
"""
import contextlib
import csv
import gzip
import io
import itertools
import json
import logging
import os

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from bedesign import ENGINE_VERSION

from .cachekey import cache_key, manifest_key, result_key
from .config import Settings
from .jobs import JobPool, PoolBusy
from .params import ValidationError, parse_designs_query
from .ratelimit import TokenBucket, client_ip
from .references import References
from .storage import LocalStorage

log = logging.getLogger(__name__)

TEMPLATES = Jinja2Templates(directory=os.path.join(os.path.dirname(__file__), 'templates'))

# How often the running page re-checks. Design doc 3 says 3 seconds.
POLL_SECONDS = 3
# Slice 1 renders a fixed slice of the table; slice 2 adds paging and filters.
PREVIEW_ROWS = 200

# No script-src at all: the running page re-checks with a meta refresh rather than
# JavaScript, so slice 1 needs no script origin. Slice 2 relaxes this to 'self' when
# it vendors HTMX for the table fragments.
CSP = ("default-src 'self'; script-src 'none'; style-src 'self'; img-src 'self' data:; "
       "connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'")
SECURITY_HEADERS = {
    'Content-Security-Policy': CSP,
    'X-Content-Type-Options': 'nosniff',
    'Referrer-Policy': 'no-referrer',
    'X-Frame-Options': 'DENY',
    'Cross-Origin-Opener-Policy': 'same-origin',
}

ROBOTS = 'User-agent: *\nDisallow: /\n'


def scrub(value):
    """Makes a user-supplied string safe to put in a log line.

    Without this, a query string containing newlines could write whatever it liked
    into the log, including entries that look like they came from the server.
    """
    return repr(str(value)[:200])


def create_app(settings=None):
    settings = settings or Settings.from_env()

    @contextlib.asynccontextmanager
    async def lifespan(app):
        app.state.references = References(settings.refdata, settings.clinvar_db)
        app.state.storage = LocalStorage(settings.results_dir)
        app.state.pool = JobPool(app.state.references, settings)
        app.state.limiter = TokenBucket(settings.rate_burst, settings.rate_seconds)
        log.info('serving %s, clinvar %s, results in %s',
                 app.state.references.bundle, app.state.references.clinvar_version,
                 app.state.storage.root)
        try:
            yield
        finally:
            app.state.pool.shutdown()

    app = FastAPI(title='Base editor guide designs', docs_url=None, redoc_url=None,
                  lifespan=lifespan)
    app.state.settings = settings
    app.mount('/static', StaticFiles(
        directory=os.path.join(os.path.dirname(__file__), 'static')), name='static')

    @app.middleware('http')
    async def security_headers(request, call_next):
        response = await call_next(request)
        for name, value in SECURITY_HEADERS.items():
            response.headers.setdefault(name, value)
        return response

    def page(request, template, status_code=200, **context):
        context.update(poll_seconds=POLL_SECONDS)
        return TEMPLATES.TemplateResponse(request=request, name=template,
                                          context=context, status_code=status_code)

    @app.get('/healthz')
    def healthz():
        return JSONResponse({'status': 'ok', 'engine_version': ENGINE_VERSION})

    @app.get('/robots.txt', response_class=PlainTextResponse)
    def robots():
        # Design doc 3 names the trade-off: a GET can start work, so crawlers are
        # turned away here and again with a noindex meta tag in base.html.
        return PlainTextResponse(ROBOTS)

    @app.get('/designs', response_class=HTMLResponse)
    def designs(request: Request):
        references = app.state.references
        try:
            transcript_id, params = parse_designs_query(
                request.query_params, references.transcript_exists)
        except ValidationError as e:
            return page(request, 'error.html', status_code=400, message=str(e))

        key = cache_key(transcript_id, params, references.release,
                        references.clinvar_version, ENGINE_VERSION)

        manifest = _manifest(app.state.storage, key)
        if manifest is not None:
            return _table(request, key, transcript_id, params, manifest)

        failed = app.state.pool.failure(key)
        if failed:
            return page(request, 'error.html', status_code=500, message=failed)

        if not app.state.pool.running(key):
            client = client_ip(request, settings.trusted_proxy)
            if not app.state.limiter.allow(client):
                return page(request, 'error.html', status_code=429,
                            message='Too many designs started from here. '
                                    'Wait a moment and reload.')
            try:
                started = app.state.pool.submit(key, transcript_id, params)
            except PoolBusy:
                response = page(request, 'busy.html', status_code=503)
                response.headers['Retry-After'] = str(POLL_SECONDS)
                return response
            if started:
                log.info('started %s for %s', key[:12], scrub(transcript_id))

        return page(request, 'running.html', transcript_id=transcript_id, key=key)

    def _table(request, key, transcript_id, params, manifest):
        header, rows = _read_designs(app.state.storage, key)
        return page(request, 'table.html', transcript_id=transcript_id, params=params,
                    key=key, header=header, rows=rows,
                    total=manifest.get('designs', len(rows)), shown=len(rows))

    return app


def _manifest(storage, key):
    """The manifest for a key, or None if this result is not cached.

    The manifest is written last, so finding it means the objects behind it are
    complete. Reading it here rather than testing for it costs the same stat and
    yields the row count, which saves parsing the whole designs file for it.
    """
    try:
        return json.loads(storage.get(manifest_key(key)))
    except KeyError:
        return None


def _read_designs(storage, key, limit=PREVIEW_ROWS):
    """The first `limit` design rows, for the slice 1 table.

    Decompressed lazily and cut off with islice: TTN is 25,662 rows, and inflating
    and parsing all of them to show 200 is work the page never uses. The row count
    comes from the manifest instead.

    Reads only a file this app wrote, at a key built from a digest.
    """
    raw = storage.get(result_key(key, 'designs'))
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
        text = io.TextIOWrapper(gz, encoding='utf-8', newline='')
        reader = csv.reader(text, delimiter='\t')
        header = next(reader, [])
        return header, list(itertools.islice(reader, limit))


app = create_app()
