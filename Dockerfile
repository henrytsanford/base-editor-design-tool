# syntax=docker/dockerfile:1
#
# The application image: source and pinned dependencies, no reference data.
# deploy/Dockerfile.bundled adds the bundle on top; deploy/README.md has the commands.
#
# Cloud Run serves linux/amd64, so a build on an arm64 machine needs
# `--platform linux/amd64` before it is pushed.

FROM python:3.13-slim

# Unbuffered so log lines reach the collector as they happen rather than when the
# block fills. No .pyc, because the layer holding them is read-only.
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# RESULTS_DIR is the only path the service itself writes to. On Cloud Run this
# filesystem is in memory, so what lands there counts against --memory.
#
# The app directory is writable too because the CLI writes its timestamped run
# folder relative to the working directory, which test/test_golden.py and
# test/test_design.py depend on.
RUN groupadd --system bedesign \
    && useradd --system --gid bedesign --create-home --home-dir /opt/bedesign \
        --shell /usr/sbin/nologin bedesign \
    && install -d -o bedesign -g bedesign /opt/bedesign/app /opt/bedesign/results

WORKDIR /opt/bedesign/app

# Dependencies are their own layer: they change far less often than the engine,
# so editing bedesign/ does not reinstall pandas.
COPY requirements.txt requirements-service.txt constraints.txt ./
RUN pip install --no-cache-dir -r requirements-service.txt -c constraints.txt

COPY --chown=bedesign:bedesign . .

# REFDATA names the parent rather than the bundle: find_bundle picks the highest
# ensembl-<release> under it and find_clinvar_db the newest clinvar-<date>.db, so
# the image does not have to be edited when the bundle is rebuilt.
ENV REFDATA=/opt/bedesign/refdata \
    RESULTS_DIR=/opt/bedesign/results

USER bedesign
EXPOSE 8000

# Ignored by Cloud Run, which has health checking of its own; this is for the
# person running the image by hand.
HEALTHCHECK --interval=30s --timeout=3s --start-period=10s CMD \
    python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:${PORT:-8000}/healthz', timeout=2)"

# Binding 0.0.0.0 is safe here because the container's network namespace is the
# boundary and the platform supplies $PORT. `exec` so that uvicorn is PID 1 and
# receives SIGTERM directly.
CMD ["sh", "-c", "exec uvicorn service.app:app --host 0.0.0.0 --port ${PORT:-8000}"]
