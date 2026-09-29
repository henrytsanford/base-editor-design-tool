"""The seven tools, as plain functions over plain dicts.

Nothing here imports MCP. The server module turns these into tool declarations and
maps their exceptions onto protocol errors; keeping the bodies free of that means the
whole surface can be called and tested directly, and that a change of transport is a
change to one file.

Every argument that service/params.py already knows how to check goes through it, via
the Query adapter. The PAM pattern, the window rule, the sg_len and intron_buffer
ranges, the preset conflict and the filter vocabularies are inherited from the web
app rather than restated, so the two front doors cannot disagree about what a valid
request is.
"""
import base64
import dataclasses
import gzip
import itertools
import re
import time

import anyio

from bedesign import ANNOTATION_COLUMNS, ENGINE_VERSION, tsv
from bedesign.engine import BE_TYPES, DEFAULT_BE_TYPE, DesignParams

from service import params as svc_params
from service.params import ValidationError
from service.references import GENE_LIMIT, TRANSCRIPT_LIMIT
from service.cachekey import RESULT_FILES, Target, cache_key
from service.jobs import PoolBusy
from service.ratelimit import client_ip
from service.results import ROWS_PER_PAGE, page_of

from . import summary
from .query import Query
from .runs import DesignFailed, ServiceBusy, filename_for

RUN_ID = re.compile(r'^[0-9a-f]{64}$')
# The nucleotide path marks introns by lower case, so case is meaningful and the
# sequence is never normalised.
SEQUENCE = re.compile(r'^[ACGTacgt]+$')
# A designed guide, at whatever lengths this engine will design. Built from the same
# range params.py checks sg_len against, so widening that cannot make a guide the
# server just designed unrecognisable here.
GUIDE = re.compile(r'^[ACGT]{%d,%d}$' % svc_params.SG_LEN_RANGE)
# Long enough for any transcript someone would paste, short enough that a bad paste
# is refused before the engine spends a core on it.
MAX_SEQUENCE = 100000
# A text export past this is refused rather than truncated: half a TSV that looks
# whole is worse than an error naming the other encoding.
DEFAULT_MAX_BYTES = 5000000
# The most a caller may raise it to.
MAX_EXPORT_BYTES = 20000000


@dataclasses.dataclass
class AppContext:
    """Everything a tool call needs, built once by the server's lifespan."""
    settings: object
    references: object
    runs: object
    # The web app's, when this server is mounted on it. None on stdio, which is the
    # signal to design in-process and to charge nobody.
    pool: object = None
    limiter: object = None


def _run_id(value):
    if not RUN_ID.match(value or ''):
        raise ValidationError('run_id must be the handle design_guides returned.')
    return value


# ---- reference lookups -------------------------------------------------------

def resolve_gene(ctx, symbol):
    """A gene symbol, resolved against the bundle."""
    resolved, _ = svc_params.parse_genes_query(Query({'q': symbol}))
    gene, matches = ctx.references.resolve_gene(resolved)
    return {'query': resolved,
            'gene': gene,
            'matches': matches,
            'truncated': len(matches) >= GENE_LIMIT}


def list_transcripts(ctx, gene, limit=TRANSCRIPT_LIMIT):
    """Every transcript of a gene, best choice first."""
    resolved, _ = svc_params.parse_genes_query(Query({'q': gene}))
    if not resolved:
        raise ValidationError('A gene symbol is required, e.g. MAP2K1.')
    limit = max(1, min(int(limit), TRANSCRIPT_LIMIT))
    found = ctx.references.transcripts_for_gene(resolved, limit)
    for i, transcript in enumerate(found):
        transcript['recommended'] = i == 0
    return {'gene': resolved, 'count': len(found),
            'truncated': len(found) >= limit, 'transcripts': found}


def list_editors(ctx):
    """The base editors this engine knows, and the bounds on a custom one.

    Built from the engine's own preset table, so this can never offer an editor the
    validator would refuse.
    """
    editors = []
    for name in BE_TYPES:
        resolved = DesignParams.from_preset(name)
        editors.append({'name': name, 'pam': resolved.pam, 'window': resolved.window,
                        'sg_len': resolved.sg_len, 'edit': resolved.edit})
    return {
        'default': DEFAULT_BE_TYPE,
        'editors': editors,
        'edits': list(svc_params.EDITS),
        'engine_version': ENGINE_VERSION,
        'limits': {
            'sg_len': list(svc_params.SG_LEN_RANGE),
            'intron_buffer': list(svc_params.INTRON_BUFFER_RANGE),
            'pam': '2-8 IUPAC nucleotide codes, e.g. NGG',
            'window': 'start-end, with 1 <= start <= end <= sg_len',
        },
    }


# ---- designing ---------------------------------------------------------------

async def design_guides(ctx, request=None, transcript_id=None, sequence=None,
                        sequence_name=None, preset=None, pam=None, window=None,
                        sg_len=None, edit=None, intron_buffer=None, filter_gc=None):
    """Designs every guide over a transcript or a pasted sequence.

    Returns the shape of the result and a handle to it, never the guides: a run is
    routinely thousands of rows, and which of them matter is the next question rather
    than this one.
    """
    editor = {'preset': preset, 'pam': pam, 'window': window, 'sg_len': sg_len,
              'edit': edit, 'intron_buffer': intron_buffer, 'filter_gc': filter_gc}
    if (transcript_id is None) == (sequence is None):
        raise ValidationError(
            'Give either transcript_id or sequence, not both and not neither.')

    if transcript_id is not None:
        resolved, design_params = svc_params.parse_designs_query(
            Query(dict(editor, transcript=transcript_id)),
            ctx.references.transcript_exists)
        target = Target(kind='transcript', transcript_id=resolved)
    else:
        design_params = svc_params.parse_editor_params(Query(editor))
        target = Target(kind='sequence', name=_sequence_name(sequence_name),
                        sequence=_sequence(sequence))

    if ctx.pool is None:
        # stdio: no pool, so the design runs in this process -- on a worker thread,
        # which is where the SDK already ran this whole function before it had to
        # wait on anything.
        run_id, manifest = await anyio.to_thread.run_sync(
            ctx.runs.design, target, design_params)
    else:
        run_id = cache_key(target.key, design_params, ctx.references.release,
                           ctx.references.clinvar_version, ENGINE_VERSION)
        manifest = await _pooled(ctx, run_id, target, design_params, request)
    table = None
    if manifest['designs'] <= ctx.settings.table_max_rows:
        # Off the event loop: a 50,000-row frame is ~65 MB of pandas work, and this
        # is the one tool that does not already run on a thread.
        table = await anyio.to_thread.run_sync(ctx.runs.table, run_id)
    return summary.build(manifest, table)


# How often a waiting design_guides looks for the manifest the worker writes last.
# Designs run from well under a second to the job timeout, so this is noise.
POLL_SECONDS = 0.25


async def _pooled(ctx, run_id, target, params, request):
    """The design, run in the web app's process pool and waited for here.

    The same pool /designs uses, which is the point: a design started in a browser and
    the same design asked for here are one job writing one result, under the cache key
    service/cachekey.py computes for both. Nothing comes back through the pool -- the
    worker writes the files and then the manifest, so the manifest's presence is what
    says the run is complete.
    """
    manifest = ctx.runs.lookup(run_id)
    if manifest is not None:
        return manifest                 # cached: free, and never charged a token
    if ctx.pool.failure(run_id):
        raise DesignFailed(_failure_message(ctx, run_id))
    await _start(ctx, run_id, target, params, request)
    return await _wait(ctx, run_id)


async def _start(ctx, run_id, target, params, request):
    """Starts the job, or joins one already running for this key.

    /designs' order, for /designs' reasons: a running job is joined for free, and only
    a call that would start one is charged a rate-limit token.

    `client=None` is deliberate. The pool's one-job-per-client rule stops a browser
    polling at its rate limit from holding every worker, but keyed on an address here
    it would serialise a whole NAT, or a hosted model's egress, to one design at a
    time -- with no page to explain the wait. What bounds this door instead is
    PoolBusy and the token bucket.
    """
    if ctx.pool.running(run_id):
        return
    client = (client_ip(request, ctx.settings.trusted_proxy_hops)
              if request is not None else None)
    if client is not None and ctx.limiter is not None and not ctx.limiter.allow(client):
        raise ServiceBusy(
            'Too many designs started from this address. Wait about %d seconds and '
            'call design_guides again with the same arguments. Reading a run you '
            'already have, with query_guides or export_run, is not limited.'
            % ctx.settings.rate_seconds)
    deadline = time.monotonic() + ctx.settings.mcp_wait
    while not ctx.pool.running(run_id):
        try:
            ctx.pool.submit(run_id, target, params)
            return
        except PoolBusy:
            # There is no queue, by design: /designs answers a busy pool with a page
            # that retries. Waiting here is better than spending the model's turn,
            # but only for as long as the token this call was already charged is
            # worth.
            if time.monotonic() >= deadline:
                raise ServiceBusy(
                    'Every design worker is busy and this design was not started. '
                    'Nothing was lost: call design_guides again with the same '
                    'arguments in about a minute.')
            await anyio.sleep(POLL_SECONDS)


async def _wait(ctx, run_id):
    """Waits for the manifest, which the worker writes last.

    Polled rather than awaited on a future: submit() reports whether this call started
    the job, not which future runs it, and a job /designs started a moment earlier may
    have no future left to hand over.

    `running` is sampled before the manifest is read, because the worker writes the
    manifest, returns, and only then is dropped from the running table -- so a job
    that finishes between the two reads is seen as finished on this pass rather than
    declared gone.
    """
    deadline = time.monotonic() + ctx.settings.job_timeout + ctx.settings.mcp_wait
    while True:
        running = ctx.pool.running(run_id)
        manifest = ctx.runs.lookup(run_id)
        if manifest is not None:
            return manifest
        if ctx.pool.failure(run_id):
            raise DesignFailed(_failure_message(ctx, run_id))
        if not running:
            raise DesignFailed(
                'The result was cleared before it could be read. Call design_guides '
                'again with the same arguments.')
        if time.monotonic() >= deadline:
            raise ServiceBusy(
                'This design is still running. Call design_guides again with exactly '
                'the same arguments -- it will pick up the finished result rather '
                'than starting over.')
        await anyio.sleep(POLL_SECONDS)


def _failure_message(ctx, run_id):
    """The pool's failure, said in terms a model can act on.

    The pool's own message ends 'Reload to try again', which is advice for a browser.
    """
    if ctx.pool.failure_reason(run_id) == 'timeout':
        return ('This design ran longer than %d seconds and was stopped. Narrow it: '
                'a more specific PAM, one deaminase (edit="C-T" or edit="A-G") '
                'rather than "all", or a shorter transcript.'
                % ctx.settings.job_timeout)
    return ('This design failed on the server and no result was stored. Call '
            'design_guides again; if it fails a second time the transcript or the '
            'parameter combination is at fault rather than the service.')


def _sequence(value):
    if len(value) > MAX_SEQUENCE:
        raise ValidationError('sequence must be at most %d bases.' % MAX_SEQUENCE)
    if not SEQUENCE.match(value):
        raise ValidationError(
            'sequence must be nucleotides: A, C, G or T. Lower case marks intronic '
            'bases, as it does in the FASTA input.')
    return value


def _sequence_name(value):
    if not value:
        raise ValidationError('sequence_name is required alongside sequence.')
    if not svc_params.GENE_QUERY.match(value):
        raise ValidationError(
            'sequence_name is letters, digits, dot, dash or underscore, e.g. GFP.')
    return value


# ---- reading a run back ------------------------------------------------------

def query_guides(ctx, run_id, page=1, consequence=None, mutation=None,
                 significance=None, deaminase=None, strand=None, hide_bsmbi=False,
                 hide_4t=False, exon=None, sub=None, sort=None, dir='asc',
                 columns=None):
    """One page of designed guides, filtered and sorted."""
    run_id = _run_id(run_id)
    view = svc_params.parse_view_query(Query({
        'page': page, 'consequence': consequence, 'mutation': mutation,
        'significance': significance, 'deaminase': deaminase, 'strand': strand,
        'hide_bsmbi': hide_bsmbi, 'hide_4t': hide_4t, 'exon': exon, 'sub': sub,
        'sort': sort, 'dir': dir}))
    table = ctx.runs.table(run_id)
    selected = table.select(view)

    names = list(table.columns)
    rows = selected.rows
    if columns:
        unknown = [name for name in columns if name not in names]
        if unknown:
            raise ValidationError('Not a column of the table: %s.'
                                  % ', '.join(unknown))
        keep = [names.index(name) for name in columns]
        names = list(columns)
        rows = [[row[i] for i in keep] for row in rows]

    return {'run_id': run_id, 'columns': names, 'rows': rows,
            'matched': selected.matched, 'page': selected.page,
            'pages': selected.pages, 'rows_per_page': ROWS_PER_PAGE,
            # Derived from the view the same way TableView.filtered is, so a filter
            # added there is echoed back here without being listed twice.
            'filters': {field.name: getattr(view, field.name)
                        for field in dataclasses.fields(view)
                        if field.name not in ('sort', 'dir', 'page')
                        and getattr(view, field.name)},
            'sort': view.sort, 'dir': view.dir}


def get_clinvar_annotations(ctx, run_id, page=1, sgrna=None, matched_only=False):
    """One page of the ClinVar annotation file.

    Rows are ragged on purpose: a guide edit that matched a ClinVar SNP carries all 28
    fields, and one that matched nothing stops after 'Mutation category'. They are not
    padded, because a blank 'SNP clinical significance' and no SNP at all are
    different answers.
    """
    run_id = _run_id(run_id)
    # Through the same validator the view parameters use, so `page` is bounded and
    # refused in one place rather than two.
    page = svc_params.parse_view_query(Query({'page': page})).page
    manifest = ctx.runs.manifest(run_id)
    wanted = None
    if sgrna:
        wanted = sgrna.upper()
        if not GUIDE.match(wanted):
            raise ValidationError('sgrna must be a guide sequence, %d-%d bases.'
                                  % svc_params.SG_LEN_RANGE)

    keep = _annotation_filter(wanted, matched_only)
    with ctx.runs.annotation_reader(manifest) as reader:
        if keep is None:
            # Nothing to count: the manifest already knows how many rows there are,
            # so the page can be sliced straight out of the stream.
            matched = manifest['annotations']
            page, pages, start = page_of(matched, page)
            rows = list(itertools.islice(reader, start, start + ROWS_PER_PAGE))
        else:
            # A filter has to see every row, but only the page is kept.
            matched, rows = 0, []
            for row in reader:
                if not keep(row):
                    continue
                matched += 1
                rows.append(row)
            page, pages, start = page_of(matched, page)
            rows = rows[start:start + ROWS_PER_PAGE]

    return {'run_id': run_id, 'columns': list(ANNOTATION_COLUMNS),
            'rows': rows, 'ragged': True,
            'matched': matched, 'page': page, 'pages': pages,
            'rows_per_page': ROWS_PER_PAGE}


def _annotation_filter(sgrna, matched_only):
    """What to keep, or None when every row qualifies."""
    if not sgrna and not matched_only:
        return None
    full = len(ANNOTATION_COLUMNS)

    def keep(row):
        if sgrna and (not row or row[0] != sgrna):
            return False
        return not matched_only or len(row) == full

    return keep


def export_run(ctx, run_id, file='designs', encoding='text',
               max_bytes=DEFAULT_MAX_BYTES):
    """A whole output file, exactly as the CLI would have written it.

    'text' is the file's own text. 'base64-gzip' is the stored gzip, for a result too
    large to read inline. Nothing is ever truncated: a partial TSV that looks whole is
    worse than being told to ask for the other encoding.
    """
    run_id = _run_id(run_id)
    if file not in RESULT_FILES:
        raise ValidationError('file must be one of %s.'
                              % ', '.join(sorted(RESULT_FILES)))
    manifest = ctx.runs.manifest(run_id)
    # Size, lines and digest were measured when the file was written, so an oversize
    # export is refused without decompressing the file it is refusing.
    stats = ctx.runs.stats(manifest, file)
    payload = dict(stats, run_id=run_id, file=file, encoding=encoding,
                   filename=filename_for(file, manifest['label']))
    if encoding == 'text':
        if stats['bytes'] > max_bytes:
            raise ValidationError(
                '%s is %d bytes, over the %d asked for. Raise max_bytes, or ask for '
                "encoding='base64-gzip'." % (payload['filename'], stats['bytes'],
                                             max_bytes))
        payload['data'] = tsv.as_text(gzip.decompress(ctx.runs.stored(manifest, file)))
    else:
        # The stored bytes, handed over without being decompressed to describe them.
        payload['data'] = base64.b64encode(
            ctx.runs.stored(manifest, file)).decode('ascii')
    return payload
