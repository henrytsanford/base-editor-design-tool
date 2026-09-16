"""Settings, read from the environment.

Everything here is explicit and local. There are no secrets, and no default that
points outside the repository -- the S3 settings the design doc lists (BUCKET,
ENSEMBL_RELEASE, CLINVAR_VERSION) arrive with S3Storage at M2.

Numeric settings are clamped rather than trusted. A MAX_JOBS of 10000 in a stray
environment variable would otherwise fork until the machine died.
"""
import os
from dataclasses import dataclass

# Upper bounds exist so a typo in the environment cannot exhaust the host.
MAX_JOBS_LIMIT = 16
JOB_TIMEOUT_LIMIT = 3600


def _bounded_int(env, name, default, low, high):
    raw = env.get(name)
    if raw is None or raw == '':
        return default
    try:
        value = int(raw)
    except ValueError:
        raise ValueError('%s must be an integer, got %r' % (name, raw))
    if not low <= value <= high:
        raise ValueError('%s must be between %d and %d, got %d' % (name, low, high, value))
    return value


def _flag(env, name, default=False):
    raw = env.get(name)
    if raw is None or raw == '':
        return default
    return raw.strip().lower() in ('1', 'true', 'yes', 'on')


@dataclass(frozen=True)
class Settings:
    refdata: str = 'refdata'
    clinvar_db: str = ''
    results_dir: str = 'results'
    max_jobs: int = 2
    # Wall-clock cap per design job. TTN, the worst gene in the reference set,
    # takes 43 s (design doc 4.3), so this leaves generous headroom while still
    # bounding a pathological request.
    job_timeout: int = 300
    # How long a failed key is remembered before a reload retries it.
    failure_ttl: int = 300
    # Token bucket on the path that starts a job: `rate_burst` immediately, then
    # one more every `rate_seconds`.
    rate_burst: int = 5
    rate_seconds: int = 30
    # Off by default and deliberately so: see trust_client_ip below.
    trusted_proxy: bool = False

    @classmethod
    def from_env(cls, env=None):
        env = os.environ if env is None else env
        return cls(
            refdata=env.get('REFDATA') or 'refdata',
            clinvar_db=env.get('CLINVAR_DB') or '',
            results_dir=env.get('RESULTS_DIR') or 'results',
            max_jobs=_bounded_int(env, 'MAX_JOBS', 2, 1, MAX_JOBS_LIMIT),
            job_timeout=_bounded_int(env, 'JOB_TIMEOUT', 300, 1, JOB_TIMEOUT_LIMIT),
            failure_ttl=_bounded_int(env, 'FAILURE_TTL', 300, 0, 86400),
            rate_burst=_bounded_int(env, 'RATE_BURST', 5, 1, 1000),
            rate_seconds=_bounded_int(env, 'RATE_SECONDS', 30, 1, 3600),
            trusted_proxy=_flag(env, 'TRUSTED_PROXY'),
        )
