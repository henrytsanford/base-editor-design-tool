# Base editor design tool
<p>This code designs every possible guide in the CDS region, 30 nucleotides into the intronic and UTR regions for 
user-defined transcripts. It then annotates the possible edits for each guide. Separate files annotating ClinVar SNPs are also generated.</p>
<b>Author</b>: Mudra Hegde, Ruth Hanna <br/>
<b>Email</b>: mhegde@broadinstitute.org, rhanna@broadinstitute.org <br/>
<b>Version: 3.0 </b> 

<b>Inputs</b>
1. <b>Input File</b>:.txt file with list of Ensembl transcript IDs in the first column and gene symbols in the second column OR 
FASTA file with nucleotide sequence
2. <b>ClinVar database</b>: ClinVar SNPs, as a database built by <code>tools/build_clinvar.py</code>; Default: the newest <code>clinvar-&lt;date&gt;.db</code> under <code>--refdata</code>. Pass <code>--no-clinvar</code> to skip the annotation.
3. <b>Input type</b>: Indicate whether the file contains a list of transcripts or a nucleotide sequence.
4. <b>Base editor type</b>: Indicate the type of base editor (Rees et al.,2018) for which designs are required. This choice dictates the choice of PAM, sgRNA length, editing window and type of edit.
5. <b>PAM</b>: PAM preference if BE type has not been selected; Default: NGG.
6. <b>Edit window</b>: Editing window relative to nucleotide position in sgRNA, if BE type has not been selected; Default: 4-8.
7. <b>sgRNA length</b>: Length of sgRNA excluding PAM sequence, if BE type has not been selected; Default:20.
8. <b>Edit</b>: Type of edit made by base editor, if BE type has not been selected; Default: all, annotates for both C->T and A->G edit.
9. <b>Intron buffer</b>: Number of bp into the intron to consider for guide design.
10. <b>Filter GC</b>: Whether to filter out edits in a GC motif.
11. <b>Output name</b>: Name for output folder.
12. <b>Source</b>: Where transcript data comes from: <code>rest</code> (the Ensembl REST API, the default) or <code>local</code> (a reference bundle built by <code>tools/build_reference.py</code>). See "Working offline" below.
13. <b>Refdata</b>: Directory holding the local reference bundle and the ClinVar database; Default: refdata.


## Requirements

Python 3.9 or newer (tested on 3.13).

## Setup

    git clone https://github.com/mhegde/base-editor-design-tool.git
    cd base-editor-design-tool
    pip install -r requirements.txt

Build the ClinVar database, which the tool uses to annotate each edit with the SNPs it
would create:

    python tools/build_clinvar.py       # ~440 MB downloaded; the build itself takes ~30 s

This writes `refdata/clinvar-<date>.db`, which design runs query one gene at a time.
It populates the `Clinical significance` column and the `clinvar_annotations_*.txt`
file. NCBI publishes a new `variant_summary.txt` weekly; re-run to refresh.

## Working offline (recommended)

By default the tool fetches transcript data from `rest.ensembl.org`, five calls per
transcript. That API is frequently overloaded, and when it is down a run cannot produce
any output. You can instead build a local copy of the reference data once and design
against it with no network at all:

    python tools/build_reference.py            # ~30 minutes, ~1 GB downloaded

If a proxy mishandles the build's parallel downloads, add `--streams 1`. The finished
bundle is about 2 GB.

Then add `--source local` to any design command:

    python base_editing_guide_designs.py \
        --input-file my_transcripts.txt \
        --input-type tid \
        --source local \
        --output-name my_designs

Besides removing the outage risk this is considerably faster, and it makes results
reproducible: designs are pinned to one Ensembl release rather than to whatever the API
happened to serve that day. The release used is recorded in each run's `README.txt`.

Rebuild when you want a newer Ensembl release (`--release`); old bundles can be kept
alongside so past results stay reproducible.

## Web service (in progress)

There is also a small web app that serves the same designs over HTTP. It needs a local
reference bundle and a ClinVar database, as above, and its own dependencies:

    pip install -r requirements-service.txt
    python -m uvicorn service.app:app --port 8000

Then open a URL that spells out the request:

    http://localhost:8000/designs?transcript=ENST00000307102&preset=ABE7.10
    http://localhost:8000/designs?transcript=ENST00000307102&pam=NGG&window=4-8&sg_len=20&edit=all

The first request starts the design and shows a page that re-checks every few seconds;
the result is then cached under `results/`, keyed by a hash of the parameters, so the
same link is served from disk afterwards. Settings come from the environment:
`REFDATA`, `CLINVAR_DB`, `RESULTS_DIR`, `MAX_JOBS`, `JOB_TIMEOUT`.

The gene search page, table filtering and TSV download are not built yet.

## Example

    python base_editing_guide_designs.py \
        --input-file Sample_data/GFP.fasta \
        --input-type nuc \
        --edit C-T \
        --output-name GFP

Results are written to `GFP_<timestamp>/`.

## Troubleshooting

**`Ensembl REST API is unavailable`** — `rest.ensembl.org` is down or overloaded;
requests have already been retried seven times. Re-run later, or use `--source local`.
To check the API, request a real record: <https://rest.ensembl.org/info/ping> answers
even while data endpoints fail.

**`Transcript '...' not found in Ensembl`** — check the ID at <https://www.ensembl.org>
(version suffixes such as `.13` are stripped automatically). With `--source local`, the
ID may also be retired in, or newer than, the bundle's Ensembl release.

## Developers only

### Tests

    pytest test/

Tests that query the Ensembl REST API run only with `ENSEMBL_TESTS=1`, and tests that
need the local bundle skip until it is built. `test/test_local_source.py` checks that the
bundle and the REST API agree on a panel of transcripts covering strand, exon count,
missing UTRs and non-coding transcripts.

`test/test_golden.py` compares whole output files against frozen fixtures in
`test/golden/`, over a panel chosen to reach every structural case: both strands, a
single-exon transcript, one with no UTR, a non-coding one, and FASTA input. After an
intentional change to the output, rewrite them and review the diff:

    pytest test/test_golden.py --regenerate-golden

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
`ANNOTATION_COLUMNS` — the same content the CLI writes to its three files. Pass
`clinvar=None` to skip annotation, `DesignParams.from_preset('ABE7.10')` to use a base
editor by name (it raises `ValueError` for an unknown one), and `design_sequence(name,
sequence, params)` for raw nucleotide input. Nothing is held in module state, so a
caller can keep several sources open at once.
