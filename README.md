# Base editor design tool

Designs every possible CRISPR base-editor guide across a transcript's coding sequence, plus 30 bp
into introns and UTRs, annotates what each edit does to the protein, and cross-references ClinVar
for known variants a guide would recreate. Input is Ensembl transcript IDs or a FASTA sequence.

## About this fork

This is a fork of [mhegde/base-editor-design-tool](https://github.com/mhegde/base-editor-design-tool)
by Mudra Hegde and Ruth Hanna (Broad Institute). The guide-design algorithm is theirs; nearly
everything around it is new. This fork adds:

- **A hosted web app** at <https://bedesigner.web.app> — nothing to install.
- **An MCP server**, so a chat assistant such as Claude can design and query guides.
- **Offline reference data**: a local, release-pinned Ensembl bundle, so runs don't depend on the
  Ensembl REST API and are reproducible.
- **A local ClinVar database** built from NCBI's weekly `variant_summary.txt`.
- **An importable engine** (`bedesign/`) and golden-file tests that pin output byte for byte.
## Web app

Open <https://bedesigner.web.app>, search for a gene, pick a transcript and a base editor, and you
get a filterable table with TSV downloads. Every view is spelled out in its URL, so a result is a
link you can send to someone:

    https://bedesigner.web.app/designs?transcript=ENST00000307102&preset=ABE7.10
    https://bedesigner.web.app/designs?transcript=ENST00000307102&pam=NGG&window=4-8&sg_len=20&edit=all
    https://bedesigner.web.app/designs?transcript=ENST00000307102&preset=BE4max&mutation=Nonsense&sort=%23+edits&dir=desc

The first request starts the design and the page re-checks every few seconds; after that the
result is cached, so the same link loads immediately. Downloads are always the complete file, not
the filtered view.

## Chat assistant (MCP)

The same engine is available over the [Model Context Protocol](https://modelcontextprotocol.io),
so an assistant can plan a base-editing experiment conversationally: resolve a gene, pick the MANE
transcript, choose a deaminase, design, then ask for the guides that matter. Nothing needs to be
installed locally:

    claude mcp add --transport http bedesign https://bedesign-ptj3oo2b3a-uc.a.run.app/api/mcp

Use this Cloud Run address rather than `bedesigner.web.app`: Firebase Hosting cuts requests off at
60 s, less than a design is allowed to take.

Seven tools: `resolve_gene`, `list_transcripts` and `list_editors` to settle what to design;
`design_guides` to run it; `query_guides`, `get_clinvar_annotations` and `export_run` to read the
result back. `design_guides` returns counts, per-exon coverage and a run handle rather than rows,
since a transcript routinely yields thousands of guides.

## Command line

Requires Python 3.9 or newer (tested on 3.13).

    git clone https://github.com/henrytsanford/base-editor-design-tool.git
    cd base-editor-design-tool
    pip install -r requirements.txt
    python tools/build_clinvar.py       # ~440 MB download, ~30 s build

`build_clinvar.py` writes `refdata/clinvar-<date>.db`, which fills the `Clinical significance`
column and the `clinvar_annotations_*.txt` file. NCBI updates its data weekly; re-run to refresh.

Example:

    python base_editing_guide_designs.py \
        --input-file Sample_data/GFP.fasta \
        --input-type nuc \
        --edit C-T \
        --output-name GFP

Results are written to `GFP_<timestamp>/`.

### Options

- `--input-file` -- a `.txt` of Ensembl transcript IDs (first column) and gene symbols (second),
  or a FASTA file.
- `--input-type` -- `tid` for transcript IDs, `nuc` for a nucleotide sequence.
- `--be-type` -- a base editor by name (the Rees et al. 2018 panel plus `ABE8e-SpRY`). Sets PAM,
  sgRNA length, editing window and edit type. Default: `ABE8e-SpRY`.
- `--pam`, `--edit-window`, `--sg-len`, `--edit` -- set those individually instead of
  `--be-type`. Defaults: `NNN`, `4-8`, `20`, `A-G`. `--edit all` annotates both C->T and A->G.
- `--intron-buffer` -- bp into each intron to tile. Default: 30.
- `--filter-gc` -- filter out edits in a GC motif. Default: `False`.
- `--output-name` -- name for the output folder.
- `--source` -- `rest` (Ensembl REST API, default) or `local` (a reference bundle; see below).
- `--refdata` -- directory holding the bundle and ClinVar database. Default: `refdata`.
- `--clinvar-db` / `--no-clinvar` -- use a specific ClinVar database, or skip annotation.

### Working offline (recommended)

By default transcript data comes from `rest.ensembl.org`, which is frequently overloaded; when it
is down a run produces nothing. Build a local copy once instead:

    python tools/build_reference.py     # ~30 minutes, ~1 GB download, ~2 GB on disk

(add `--streams 1` if a proxy mishandles the parallel downloads), then pass `--source local` to any
design command. This is faster and reproducible: designs are pinned to one Ensembl release,
recorded in each run's `README.txt`. Rebuild with `--release` for a newer one; old bundles can be
kept alongside.

### Troubleshooting

**`Ensembl REST API is unavailable`** — `rest.ensembl.org` is down or overloaded; requests have
already been retried seven times. Re-run later, or use `--source local`. To check the API, request
a real record: <https://rest.ensembl.org/info/ping> answers even while data endpoints fail.

**`Transcript '...' not found in Ensembl`** — check the ID at <https://www.ensembl.org> (version
suffixes such as `.13` are stripped automatically). With `--source local`, the ID may also be
retired in, or newer than, the bundle's Ensembl release.

## Running the web service yourself

The app needs a local reference bundle and a ClinVar database, as above, plus its own
dependencies:

    pip install -r requirements-service.txt
    python -m uvicorn service.app:app --port 8000

There is no JavaScript: the search page is a plain form and the table's controls are links, which
is why the app can serve `script-src 'none'`.

Settings come from the environment, and `service/config.py` is the schema that reads them: every
name, its default, and the reasoning behind it. `MAX_JOBS` and `JOB_TIMEOUT` are the two you are
most likely to set; `TRUSTED_PROXY_HOPS` is the one worth reading before you set it, since it
decides which `X-Forwarded-For` entry, if any, the rate limiter believes.

The service also ships as a container image with the reference bundle baked in; `deploy/README.md`
covers building, running and deploying it. The image pins its dependencies with `constraints.txt`,
because `test/golden/` compares output byte for byte while `requirements.txt` only states floors.
The hosted copy is public with no login: the app's own limits on how much work each visitor can
start take the place of access control.

One constraint worth knowing: the per-job wall-clock cap uses `signal.SIGALRM`, so the service
needs a Unix host. On Windows, run the container rather than Python directly — native Windows
Python has no `SIGALRM` and the cap silently cannot be enforced.

## Developers

### Tests

    pytest test/

Tests that query the Ensembl REST API run only with `ENSEMBL_TESTS=1`, and tests that need the
local bundle skip until it is built. `test/test_local_source.py` checks that the bundle and the
REST API agree on a panel of transcripts covering strand, exon count, missing UTRs and non-coding
transcripts.

`test/test_golden.py` compares whole output files against frozen fixtures in `test/golden/`, over
a panel chosen to reach every structural case: both strands, a single-exon transcript, one with no
UTR, a non-coding one, and FASTA input. After an intentional change to the output, rewrite them
and review the diff:

    pytest test/test_golden.py --regenerate-golden

### Running the MCP server locally

Against your own bundle, over stdio:

    pip install -r requirements-mcp.txt
    python -m mcp_server              # speaks MCP over stdin/stdout

`REFDATA` must point at a bundle. Register it with your client by absolute path — a stdio server
is launched in the client's working directory, not this one:

    claude mcp add bedesign-guides --scope local \
      --env PYTHONPATH=$PWD --env REFDATA=$PWD/refdata --env RESULTS_DIR=$PWD/results \
      -- $(which python) -m mcp_server

Results are content-addressed under `RESULTS_DIR`, so asking the same question twice designs once.
Mounted on the web service, `design_guides` runs in the same process pool the browser path uses,
under the same cache key — a design started in either place is one job writing one result. Because
that pool is bounded, a design can be turned away while every worker is busy; the tool says so and
says to call again, which costs nothing.

`export_run` returns a file byte for byte as the CLI writes it, and that is checked rather than
asserted: `test/test_mcp_golden.py` drives the tools over the same six cases as the golden suite
and compares against the same fixtures, and `test/test_mcp_eval.py` replays realistic design tasks
and checks that every row handed back appears in the frozen output. The conversational layer is a
way of reading the engine, not a second implementation of it.

The MCP stack is optional: without `requirements-mcp.txt` installed, the HTTP mount is simply
absent and `pytest test/ -q` skips the MCP tests rather than failing.

### Using the design engine from Python

The engine is importable, so it can be used without the command line:

```python
from bedesign import design_transcript, DesignParams
from bedesign.transcript_source import (ClinVarSource, LocalSource,
                                        find_bundle, find_clinvar_db)

source = LocalSource(find_bundle('refdata'))
clinvar = ClinVarSource(find_clinvar_db(None, 'refdata'))

designs, errors, annotations = design_transcript(
    source, clinvar, 'ENST00000307102', DesignParams(edit='all'))
```

The three returned lists are rows under `DESIGN_COLUMNS`, `ERROR_COLUMNS` and
`ANNOTATION_COLUMNS` — the same content the CLI writes to its three files. Pass `clinvar=None` to
skip annotation, `DesignParams.from_preset('ABE7.10')` to use a base editor by name (it raises
`ValueError` for an unknown one), and `design_sequence(name, sequence, params)` for raw nucleotide
input. Nothing is held in module state, so a caller can keep several sources open at once.
