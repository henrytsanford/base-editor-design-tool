# Deploying the service

There are two images: `Dockerfile` at the repository root holds the code and its pinned
dependencies, and `deploy/Dockerfile.bundled` adds `refdata/` on top. The bundled one is
what gets deployed.

The two images are separate because they have opposite requirements: the first has to
build in seconds without 2.5 GB in its context so the tests have something small to run
in, and the second has to carry the reference data it will answer from.

## Building

    docker build -t bedesign:thin .
    docker run --rm bedesign:thin pytest test/ -q

    TAG=ensembl-116-clinvar-2026-09-09
    docker build -f deploy/Dockerfile.bundled \
        --build-arg ENSEMBL_RELEASE=116 --build-arg CLINVAR_DATE=2026-09-09 \
        -t bedesign:$TAG refdata
    docker run --rm bedesign:$TAG pytest test/ -q

The second test run is the one that matters: with the bundle present the `bundle`-marked
tests stop skipping, so `test/golden/` is compared against the real reference data using
the exact dependency versions the image deploys with. A golden failure here that passes on
the host means `constraints.txt` has drifted from what the goldens were generated with --
read the diff before changing either.

Note the bundled build's context is `refdata`, not `.`, and it copies only the release and
date you name, so older builds left in `refdata/` stay out of the image. The tag names the
same two because both are cache-key inputs: the image and the results it produces describe
each other, and a rollback is redeploying the previous tag. Build it where `refdata/`
already exists rather than in CI -- rebuilding the bundle is ~30 minutes plus about a
gigabyte of downloads.

## Running

    docker run --rm -p 8000:8000 \
        -e MAX_JOBS=2 -e JOB_TIMEOUT=120 \
        bedesign:$TAG

    curl -s localhost:8000/healthz
    curl -s 'localhost:8000/designs?transcript=ENST00000307102&preset=ABE7.10' | head

The first `/designs` call starts a job and returns the running page; the second, a few
seconds later, returns the table.

`service/config.py` is the schema for the settings and carries the defaults.
`TRUSTED_PROXY` must stay `false` on Cloud Run: Google's front end *appends* to a
client-supplied `X-Forwarded-For` rather than replacing it, and `client_ip` trusts the
first entry.

The container filesystem does not outlive the container, so `RESULTS_DIR` needs a volume
if the result cache should survive a restart. Nothing is lost without one -- a result is
recomputable from its key -- but it is recomputed. Note that nothing expires that
directory: with a volume it grows until you remove results yourself, and without one it
grows in memory on Cloud Run, counting against `--memory` for as long as the instance lives.

## Deploying a change

Rebuild both images as above and redeploy the new tag; the data and the code move
together. Restarting drops in-flight jobs, which is by design: a lost job is restarted by
the next poll. If the change touches the engine, run `pytest test/ -q` first, and if the
output changed intentionally, bump `ENGINE_VERSION` so cache keys turn over.

## Watching it

Logs go to stdout. The things worth noticing, in order: `/healthz` not answering, and the
service restarting repeatedly rather than once.

One caveat on `/healthz`: it returns a static response and does not touch the process
pool, so it reports `ok` even when every worker is dead. A service that answers `/healthz`
but 500s on `/designs` is that failure, and no platform health check can see it.
