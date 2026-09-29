"""The MCP server: tool declarations, and the reference data they run against.

Declarations only. Every body lives in tools.py, so this file is about what the model
is offered and what an error looks like to it, and a change of transport touches
nothing else.

The transport is chosen at the bottom. stdio is what a desktop client launches, but
the server object is transport-agnostic: `mcp.run('streamable-http')` serves the same
seven tools over HTTP without a line of this changing.
"""
import argparse
import logging
import sys
from contextlib import asynccontextmanager
from typing import Annotated, Any, Literal

from mcp.server import MCPServer
from mcp.server.mcpserver import Context
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations
from pydantic import Field

from service.config import Settings
from service.params import (INTRON_BUFFER_RANGE, MAX_EXON, MAX_PAGE,
                            SG_LEN_RANGE, ValidationError)
from service.references import TRANSCRIPT_LIMIT
from service.references import References
from service.results import UnknownFilterValue

from . import tools
from .runs import ResultTooLarge, RunNotFound, RunStore
from .tools import DEFAULT_MAX_BYTES, MAX_EXPORT_BYTES

log = logging.getLogger(__name__)

INSTRUCTIONS = """\
Designs CRISPR base-editor guides across an Ensembl transcript's coding sequence, its
UTRs and 30 bp into each intron, annotates what each edit does to the protein, and
cross-references ClinVar for variants a guide would recreate.

A normal session: resolve_gene to get the symbol, list_transcripts to pick one (the
first is MANE Select), list_editors to choose a deaminase, design_guides to run it,
then query_guides to read the guides that matter. design_guides never returns rows --
a run is routinely thousands of guides -- it returns counts and a run_id, and every
filter query_guides takes names a value design_guides reported.

A non-zero counts.errors means the engine stopped early on that transcript, so the
designs are a partial answer.
"""

# Everything here reads reference data and writes nothing the caller can observe
# beyond its own cached result.
READ_ONLY = ToolAnnotations(read_only_hint=True, open_world_hint=False)

# Errors a model can do something about: a bad parameter, a gene that is not there, a
# handle that has been swept, a filter value this result does not contain. They come
# back as tool errors with the message, rather than crashing the server.
RECOVERABLE = (ValidationError, UnknownFilterValue, RunNotFound, ResultTooLarge)


def _call(fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except RECOVERABLE as e:
        raise ToolError(str(e.args[0] if e.args else e))


def build_server(settings=None):
    """The server, with its reference data opened once for the process."""
    settings = settings or Settings.from_env()

    @asynccontextmanager
    async def lifespan(_server):
        # References opens the bundle and resolves the ClinVar database, so a
        # missing or stale one fails here rather than on the first design. The
        # engine's own sources are opened per thread by the RunStore.
        references = References(settings.refdata, settings.clinvar_db)
        log.info('bundle: %s (Ensembl %s)', references.bundle, references.release)
        log.info('clinvar: %s', references.clinvar_db)
        yield tools.AppContext(settings=settings, references=references,
                               runs=RunStore(settings, references))

    mcp = MCPServer('bedesign-guides', instructions=INSTRUCTIONS,
                    version=tools.ENGINE_VERSION, lifespan=lifespan)

    @mcp.tool(title='Resolve a gene symbol', annotations=READ_ONLY)
    def resolve_gene(
        ctx: Context,
        symbol: Annotated[str, Field(max_length=32, description=(
            'A gene symbol or the start of one, e.g. MAP2K1.'))],
    ) -> dict[str, Any]:
        """Finds the gene a symbol names, in the local Ensembl bundle.

        An exact symbol resolves, and so does a prefix that can only mean one gene.
        An empty `gene` with several `matches` means the prefix was ambiguous: show
        them and ask which was meant rather than guessing.
        """
        return _call(tools.resolve_gene, ctx.request_context.lifespan_context, symbol)

    @mcp.tool(title='List a gene\'s transcripts', annotations=READ_ONLY)
    def list_transcripts(
        ctx: Context,
        gene: Annotated[str, Field(max_length=32, description='An exact gene symbol.')],
        limit: Annotated[int, Field(ge=1, le=TRANSCRIPT_LIMIT)] = TRANSCRIPT_LIMIT,
    ) -> dict[str, Any]:
        """Every transcript of a gene, best choice first.

        The first entry is flagged `recommended`: MANE Select where the gene has one,
        otherwise Ensembl canonical, with protein-coding above the retained-intron and
        NMD entries a base-editor screen rarely wants. Prefer it unless the user has
        named a transcript themselves.
        """
        return _call(tools.list_transcripts,
                     ctx.request_context.lifespan_context, gene, limit)

    @mcp.tool(title='List base editors', annotations=READ_ONLY)
    def list_editors(ctx: Context) -> dict[str, Any]:
        """The base editors this engine knows, with each one's PAM and edit window.

        `edit` is which deaminase the editor carries: 'C-T' for a cytosine base
        editor, 'A-G' for an adenine one. Designing with edit='all' runs both passes,
        which is what to do when the user has not committed to an editor yet.
        """
        return _call(tools.list_editors, ctx.request_context.lifespan_context)

    @mcp.tool(title='Design guides', annotations=READ_ONLY)
    def design_guides(
        ctx: Context,
        transcript_id: Annotated[str | None, Field(description=(
            'An Ensembl transcript, e.g. ENST00000307102. Give this or sequence.'))] = None,
        sequence: Annotated[str | None, Field(description=(
            'A raw nucleotide sequence to design over instead of a transcript. Lower '
            'case marks intronic bases; case is preserved.'))] = None,
        sequence_name: Annotated[str | None, Field(max_length=32, description=(
            'What to call the sequence in the output. Required with sequence.'))] = None,
        preset: Annotated[str | None, Field(max_length=32, description=(
            'A base editor from list_editors. Cannot be combined with pam, window, '
            'sg_len or edit.'))] = None,
        pam: Annotated[str | None, Field(max_length=8, description=(
            '2-8 IUPAC codes, e.g. NGG.'))] = None,
        window: Annotated[str | None, Field(max_length=8, description=(
            'The edit window within the guide, e.g. 4-8.'))] = None,
        sg_len: Annotated[int | None, Field(ge=SG_LEN_RANGE[0],
                                    le=SG_LEN_RANGE[1])] = None,
        edit: Annotated[Literal['C-T', 'A-G', 'all'] | None, Field(description=(
            'Which deaminase pass to run.'))] = None,
        intron_buffer: Annotated[int | None, Field(
            ge=INTRON_BUFFER_RANGE[0], le=INTRON_BUFFER_RANGE[1], description=(
            'How far into each intron to design, in bp. Allowed alongside a preset.'))] = None,
        filter_gc: Annotated[bool | None, Field(description=(
            'Drop guides with a GC motif. Allowed alongside a preset.'))] = None,
    ) -> dict[str, Any]:
        """Designs every possible guide over a transcript or a pasted sequence.

        Returns the shape of the result and a `run_id`, never the guides themselves:
        a transcript routinely yields thousands, and which of them matter is the next
        question. Read them with query_guides, or take the whole file with export_run.

        The summary reports what the result actually contains -- consequence classes,
        the ClinVar classifications present, per-exon counts, the busiest residue
        substitutions -- and every query_guides filter names one of those values.

        Watch counts.errors: a non-zero count means the engine abandoned the
        transcript part-way, so the designs are incomplete rather than final.
        """
        return _call(tools.design_guides, ctx.request_context.lifespan_context,
                     transcript_id=transcript_id, sequence=sequence,
                     sequence_name=sequence_name, preset=preset, pam=pam,
                     window=window, sg_len=sg_len, edit=edit,
                     intron_buffer=intron_buffer, filter_gc=filter_gc)

    @mcp.tool(title='Query designed guides', annotations=READ_ONLY)
    def query_guides(
        ctx: Context,
        run_id: Annotated[str, Field(description='The handle design_guides returned.')],
        page: Annotated[int, Field(ge=1, le=MAX_PAGE)] = 1,
        consequence: Annotated[
            Literal['lof', 'mis', 'sil', 'nc', 'none'] | None, Field(description=(
                "The guide's worst consequence: lof (nonsense or a broken splice "
                'site), mis (missense), sil (silent), nc (UTR, intron or flanking), '
                'none (no edit in the window).'))] = None,
        mutation: Annotated[str | None, Field(max_length=64, description=(
            "One of the engine's own categories, e.g. Nonsense or Splice-donor."))] = None,
        significance: Annotated[str | None, Field(max_length=64, description=(
            "A ClinVar classification as design_guides reported it, e.g. Pathogenic, "
            "or 'any-match' for every guide that recreates a classified variant."))] = None,
        deaminase: Annotated[Literal['C-T', 'A-G'] | None, Field(description=(
            'Keep only guides from one deaminase pass.'))] = None,
        strand: Literal['sense', 'antisense'] | None = None,
        hide_bsmbi: Annotated[bool, Field(description=(
            'Drop guides carrying a BsmBI site, which complicates cloning.'))] = False,
        hide_4t: Annotated[bool, Field(description=(
            'Drop guides with a TTTT run, which terminates Pol III transcription.'))] = False,
        exon: Annotated[int | None, Field(ge=1, le=MAX_EXON, description=(
            'Keep only guides in this exon, numbered from 1.'))] = None,
        sub: Annotated[str | None, Field(max_length=16, description=(
            'Keep only guides making this residue change, e.g. Glu-Gly.'))] = None,
        sort: Annotated[str | None, Field(max_length=64, description=(
            'A column name to sort by.'))] = None,
        dir: Literal['asc', 'desc'] = 'asc',
        columns: Annotated[list[str] | None, Field(description=(
            'Return only these columns. Omit for every column of the file.'))] = None,
    ) -> dict[str, Any]:
        """One page of designed guides, filtered and sorted. 50 rows a page.

        Every filter names a value design_guides reported for this run; one it did not
        report is refused rather than quietly returning nothing, so an empty page
        means no guide qualifies, not that the filter was misspelled.

        Rows are the engine's own text, in `columns` order.
        """
        return _call(tools.query_guides, ctx.request_context.lifespan_context,
                     run_id=run_id, page=page, consequence=consequence,
                     mutation=mutation, significance=significance,
                     deaminase=deaminase, strand=strand, hide_bsmbi=hide_bsmbi,
                     hide_4t=hide_4t, exon=exon, sub=sub, sort=sort, dir=dir,
                     columns=columns)

    @mcp.tool(title='Get ClinVar annotations', annotations=READ_ONLY)
    def get_clinvar_annotations(
        ctx: Context,
        run_id: Annotated[str, Field(description='The handle design_guides returned.')],
        page: Annotated[int, Field(ge=1, le=MAX_PAGE)] = 1,
        sgrna: Annotated[str | None, Field(max_length=24, description=(
            'Only rows for this exact guide sequence.'))] = None,
        matched_only: Annotated[bool, Field(description=(
            'Only rows where the edit matched a known ClinVar variant.'))] = False,
    ) -> dict[str, Any]:
        """Per-edit ClinVar detail for a run: which known variant each edit recreates.

        Rows are ragged by design. An edit that matched a ClinVar SNP carries all 28
        fields; one that matched nothing stops after 'Mutation category'. They are not
        padded, because a blank classification and no variant at all are different
        answers. Read a row by position against `columns`, and check its length before
        reading a SNP field.
        """
        return _call(tools.get_clinvar_annotations,
                     ctx.request_context.lifespan_context, run_id=run_id, page=page,
                     sgrna=sgrna, matched_only=matched_only)

    @mcp.tool(title='Export a run', annotations=READ_ONLY)
    def export_run(
        ctx: Context,
        run_id: Annotated[str, Field(description='The handle design_guides returned.')],
        file: Literal['designs', 'errors', 'clinvar'] = 'designs',
        encoding: Literal['text', 'base64-gzip'] = 'text',
        max_bytes: Annotated[int, Field(
            ge=1, le=MAX_EXPORT_BYTES)] = DEFAULT_MAX_BYTES,
    ) -> dict[str, Any]:
        """A whole output file, byte for byte as the command-line tool writes it.

        Use this when the user wants the data itself rather than an answer about it --
        to save, to hand to a collaborator, or to load into their own analysis. Too
        large for `text` raises rather than truncating; ask again with
        encoding='base64-gzip'.
        """
        return _call(tools.export_run, ctx.request_context.lifespan_context,
                     run_id=run_id, file=file, encoding=encoding,
                     max_bytes=max_bytes)

    return mcp


def main(argv=None):
    parser = argparse.ArgumentParser(description='Base-editor guide design over MCP.')
    parser.add_argument('--transport', default='stdio',
                        choices=['stdio', 'streamable-http', 'sse'],
                        help='stdio for a desktop client; the others serve over HTTP.')
    args = parser.parse_args(argv)
    # stderr, always: on stdio the other stream is the protocol.
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    build_server().run(args.transport)


if __name__ == '__main__':
    main()
