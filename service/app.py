"""The HTTP layer.

Design doc 3: the URL is the request. `GET /designs` validates, computes the cache
key, and then either renders the cached table, joins a running job, or starts one.
Nothing is stored about a request beyond the in-memory job table, so a restart costs
a recomputation rather than a broken link.

The same idea carries the rest of the app: a filter, a sort and a page are query
parameters, so every view of a result is a link someone can send. There is no
JavaScript anywhere -- the running page re-checks with a meta refresh, and sorting and
paging are ordinary anchors -- which is why `script-src` can stay `'none'`.

A request has two halves. Only the design half reaches the cache key, so filtering a
table re-renders a cached result instead of computing a new one.
"""
import contextlib
import dataclasses
import json
import logging
import os
from urllib.parse import urlencode

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from bedesign import DESIGN_COLUMNS, ENGINE_VERSION
from bedesign.engine import ALL_EDITS, BE_TYPES, DEFAULT_BE_TYPE, DesignParams

from .cachekey import RESULT_FILES, cache_key, manifest_key, result_key
from .config import Settings
from .jobs import JobPool, PoolBusy
from .params import (DESIGN_PARAMS, DOWNLOAD_PARAMS, EDITOR_ALL, EDITS, GENE_PARAMS,
                     INTRON_BUFFER_RANGE, SG_LEN_RANGE, STRANDS, VIEW_PARAMS,
                     ValidationError, parse_designs_query, parse_download_query,
                     parse_genes_query, parse_view_query)
from .ratelimit import TokenBucket, client_ip
from .references import GENE_LIMIT, References
from .results import ANY_MATCH, ResultCache, ResultTable, UnknownFilterValue
from .storage import LocalStorage

log = logging.getLogger(__name__)

TEMPLATES = Jinja2Templates(directory=os.path.join(os.path.dirname(__file__), 'templates'))

# How often the running page re-checks. Design doc 3 says 3 seconds.
POLL_SECONDS = 3

# The order parameters appear in a URL we build. Fixed, so the same view always has
# the same link and two people comparing URLs see the same string.
DESIGNS_ORDER = DESIGN_PARAMS + VIEW_PARAMS

# No script-src at all. Nothing in the app needs JavaScript: the running page uses a
# meta refresh, and the table's filters, sort and paging are links and a GET form.
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

# The preset dropdown, read from the engine's own table so it cannot offer an editor
# the validator would then refuse.
PRESETS = [dict(dataclasses.asdict(DesignParams.from_preset(name)), name=name)
           for name in BE_TYPES]


def scrub(value):
    """Makes a user-supplied string safe to put in a log line.

    Without this, a query string containing newlines could write whatever it liked
    into the log, including entries that look like they came from the server.
    """
    return repr(str(value)[:200])


def values_of(query):
    """The request's parameters as a plain dict.

    Safe only after validation, which is what rules out a repeated name; here the
    first value would win and quietly disagree with the value that was validated.
    """
    return {name: query[name] for name in query.keys()}


def build_url(path, values, order, **overrides):
    """A link to `path` carrying these parameters, with some changed.

    Names outside `order` are dropped rather than passed through, which is how a
    download link sheds the filters and a transcript link sheds the search box.
    Empty and false values are omitted so an unset filter leaves no trace in the URL.
    """
    merged = dict(values)
    merged.update(overrides)
    items = []
    for name in order:
        value = merged.get(name)
        if value is None or value == '' or value is False:
            continue
        items.append((name, 'true' if value is True else str(value)))
    return '%s?%s' % (path, urlencode(items))


def create_app(settings=None):
    settings = settings or Settings.from_env()

    @contextlib.asynccontextmanager
    async def lifespan(app):
        app.state.references = References(settings.refdata, settings.clinvar_db)
        app.state.storage = LocalStorage(settings.results_dir)
        app.state.pool = JobPool(app.state.references, settings)
        app.state.limiter = TokenBucket(settings.rate_burst, settings.rate_seconds)
        app.state.tables = ResultCache(settings.table_cache_rows,
                                       settings.table_cache_frames)
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

    def bad_request(request, message):
        return page(request, 'error.html', status_code=400, message=message)

    @app.get('/healthz')
    def healthz():
        return JSONResponse({'status': 'ok', 'engine_version': ENGINE_VERSION})

    @app.get('/robots.txt', response_class=PlainTextResponse)
    def robots():
        # Design doc 3 names the trade-off: a GET can start work, so crawlers are
        # turned away here and again with a noindex meta tag in base.html.
        return PlainTextResponse(ROBOTS)

    @app.get('/', response_class=HTMLResponse)
    def search(request: Request):
        """The search page. A plain GET form, so it needs no script and no session."""
        return page(request, 'search.html', presets=PRESETS, defaults=DesignParams(),
                    default_preset=DEFAULT_BE_TYPE, edits=EDITS,
                    sg_len_range=SG_LEN_RANGE,
                    intron_buffer_range=INTRON_BUFFER_RANGE)

    @app.get('/genes', response_class=HTMLResponse)
    def genes(request: Request):
        """Gene symbol to transcript.

        Renders matches for a prefix, and the transcript list once the symbol names
        exactly one gene. The editor choice rides along in the query string so the
        links out of here arrive at /designs fully specified.
        """
        try:
            symbol, params = parse_genes_query(request.query_params)
        except ValidationError as e:
            return bad_request(request, str(e))

        references = app.state.references
        gene, matches = references.resolve_gene(symbol)
        transcripts = references.transcripts_for_gene(gene) if gene else []

        values = values_of(request.query_params)
        for transcript in transcripts:
            transcript['url'] = build_url('/designs', values, DESIGN_PARAMS,
                                          transcript=transcript['transcript_id'])
        others = [{'name': name,
                   'url': build_url('/genes', values, GENE_PARAMS, q=name)}
                  for name in matches if name != gene]

        return page(request, 'genes.html', symbol=symbol, gene=gene,
                    transcripts=transcripts, others=others, params=params,
                    truncated=len(matches) >= GENE_LIMIT,
                    # So refining the search keeps the editor the user chose.
                    hidden=[(name, values[name]) for name in EDITOR_ALL
                            if values.get(name)])

    @app.get('/designs', response_class=HTMLResponse)
    def designs(request: Request):
        references = app.state.references
        try:
            transcript_id, params = parse_designs_query(
                request.query_params, references.transcript_exists)
            view = parse_view_query(request.query_params)
        except ValidationError as e:
            return bad_request(request, str(e))

        key = cache_key(transcript_id, params, references.release,
                        references.clinvar_version, ENGINE_VERSION)

        manifest = _manifest(app.state.storage, key)
        if manifest is not None:
            try:
                return _table(request, key, transcript_id, params, view, manifest)
            except UnknownFilterValue as e:
                return bad_request(
                    request, 'No guide in this result is annotated %s.' % e)

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

        return page(request, 'running.html', transcript_id=transcript_id)

    @app.get('/designs/download')
    def download(request: Request):
        """The complete result file, as the engine wrote it.

        Filters are a way of reading the table, not of cutting the file: a download is
        always the whole result, so it matches what the CLI produces for the same
        parameters. An uncached key is a 404 -- a download never starts a job, which
        would otherwise be a way around the rate limit on /designs.
        """
        references = app.state.references
        try:
            transcript_id, params, name = parse_download_query(
                request.query_params, references.transcript_exists)
        except ValidationError as e:
            return bad_request(request, str(e))

        key = cache_key(transcript_id, params, references.release,
                        references.clinvar_version, ENGINE_VERSION)
        try:
            data = app.state.storage.get(result_key(key, name))
        except KeyError:
            return page(request, 'error.html', status_code=404,
                        message='That result is not ready yet. Open the table first.')

        # transcript_id already matched ^ENST\\d{11}$ and the file name came out of a
        # fixed table, so neither can carry anything into the header.
        filename = '%s_%s' % (transcript_id, RESULT_FILES[name])
        return Response(content=data, media_type='application/gzip', headers={
            'Content-Disposition': 'attachment; filename="%s"' % filename})

    def _table(request, key, transcript_id, params, view, manifest):
        storage = app.state.storage

        def load():
            return ResultTable(storage.get(result_key(key, 'designs')),
                               DESIGN_COLUMNS)

        table = app.state.tables.get(key, load)
        result = table.select(view)
        values = values_of(request.query_params)

        def designs_url(**overrides):
            # Any change to what is shown returns to the first page: page 4 of a
            # freshly filtered table is usually not where the user wanted to land.
            overrides.setdefault('page', None)
            return build_url('/designs', values, DESIGNS_ORDER, **overrides)

        columns = []
        for column in table.columns:
            current = view.dir if view.sort == column else ''
            columns.append({
                'name': column,
                # Clicking the column already sorted ascending reverses it.
                'url': designs_url(sort=column,
                                   dir='desc' if current == 'asc' else 'asc'),
                'sorted': current,
            })

        pager = {
            'previous': designs_url(page=result.page - 1 if result.page > 2 else None)
                        if result.page > 1 else '',
            'next': designs_url(page=result.page + 1)
                    if result.page < result.pages else '',
        }

        return page(
            request, 'table.html', transcript_id=transcript_id, params=params,
            view=view, columns=columns, result=result,
            total=manifest.get('designs', table.total), facets=table.facets(),
            any_match=ANY_MATCH, edits=ALL_EDITS, strands=STRANDS, pager=pager,
            clear_url=build_url('/designs', values, DESIGN_PARAMS),
            hidden=[(name, values[name]) for name in DESIGN_PARAMS + ('sort', 'dir')
                    if values.get(name)],
            downloads=[(name, build_url('/designs/download', values, DOWNLOAD_PARAMS,
                                        file=name))
                       for name in RESULT_FILES])

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


app = create_app()
