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
`TRUSTED_PROXY_HOPS` is the number of proxies in front that append to
`X-Forwarded-For`; `client_ip` reads the entry that many places from the right, because
Google's front end *appends* to a client-supplied header rather than replacing it, so
everything to the left is whatever the caller sent. The default, 0, ignores the header,
which on Cloud Run puts every caller in one rate-limit bucket. Cloud Run alone should be
1, but measure it on a deployed revision before setting it. The old `TRUSTED_PROXY=true`
is refused at startup.

The container filesystem does not outlive the container, so `RESULTS_DIR` needs a volume
if the result cache should survive a restart. Nothing is lost without one -- a result is
recomputable from its key -- but it is recomputed. Either way it is kept under
`RESULTS_MAX_MB` (default 1024): after each job, the least recently viewed results are
deleted until the rest fit. On Cloud Run the directory is in memory, so that budget
counts against `--memory`.

## Access

The service is meant to be public on Cloud Run: anyone with the URL can use it, with no
login, and its own limits -- a per-client rate limit on starting designs, a bounded worker
pool, a per-job time cap -- stand in for access control. `--max-instances 1` is what
bounds the bill, since no amount of traffic can start a second instance.

Until those limits are fixed to work on Cloud Run, it is deployed with
`--no-allow-unauthenticated`, and only accounts granted `roles/run.invoker` can reach it.
What has to change first:

- **Identify clients correctly.** Done in code: `TRUSTED_PROXY_HOPS` reads
  `X-Forwarded-For` from the right. What is left is measuring the hop count on the
  deployed service and setting it (see "Running").
- **Price requests, not just count them.** Done: each client gets one running design
  at a time, and with ClinVar annotation now indexed by position, `JOB_TIMEOUT=60`
  covers every preset on the largest gene (TTN, 18 s at worst) while cutting off the
  widest custom requests.
- **Bound result size.** Done: requests that used to time out now finish, and the
  largest -- the default near-PAMless editor with `edit=all` on TTN -- would parse into
  a ~690 MB table. Results with more than `TABLE_MAX_ROWS` guides (50,000) are offered
  as downloads only, concurrent views of one result share a single parse, and the
  results folder is kept under `RESULTS_MAX_MB` (1024) by deleting the least recently
  viewed results after each job.
- **Survive a dead worker.** Done: a worker killed outright no longer leaves the pool
  unusable; it is replaced, and `/healthz` reports the crash (see "Watching it").

Then set a billing budget alert, and open it with:

    gcloud run services add-iam-policy-binding bedesign --region "$REGION" \
        --member=allUsers --role=roles/run.invoker

## Deploying a change

Rebuild both images as above and redeploy the new tag; the data and the code move
together. Restarting drops in-flight jobs, which is by design: a lost job is restarted by
the next poll. If the change touches the engine, run `pytest test/ -q` first, and if the
output changed intentionally, bump `ENGINE_VERSION` so cache keys turn over.

## Watching it

Logs go to stdout. The things worth noticing, in order: `/healthz` not answering, the
service restarting repeatedly rather than once, and -- once it is public -- the billing
budget alert. Any request keeps the instance warm, crawlers included, so a public URL
costs more than its real use would suggest; the budget alert is how you find out.

`/healthz` answers 503 `worker crashed` from a design worker dying until a job next
finishes normally. The service replaces the pool on its own, so a single 503 is not an
outage; one that persists means the pool keeps dying. Cloud Run only restarts the
container on it if an HTTP liveness probe on `/healthz` is configured -- give it a
failure threshold longer than one design, so it does not kill a job the replacement pool
is already running.
