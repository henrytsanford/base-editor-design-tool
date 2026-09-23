"""Settings, read from the environment.

Everything here is explicit and local. There are no secrets, and no default that
points outside the repository.

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
    # The parsed-result cache (design doc 5.3). Bounded by rows rather than entries
    # because a frame costs ~1.3 KB per row: TTN's 25,662 guides are ~34 MB, so a
    # count of frames alone would not bound memory at all.
    table_cache_rows: int = 100000
    table_cache_frames: int = 8
    # Off by default and deliberately so: see trust_client_ip below.
    trusted_proxy: bool = False

    @classmethod
    def from_env(cls, env=None):
        # Defaults are read off the fields above rather than repeated here, so
        # Settings() and Settings.from_env() cannot drift apart.
        env = os.environ if env is None else env
        return cls(
            refdata=env.get('REFDATA') or cls.refdata,
            clinvar_db=env.get('CLINVAR_DB') or cls.clinvar_db,
            results_dir=env.get('RESULTS_DIR') or cls.results_dir,
            max_jobs=_bounded_int(env, 'MAX_JOBS', cls.max_jobs, 1, MAX_JOBS_LIMIT),
            job_timeout=_bounded_int(env, 'JOB_TIMEOUT', cls.job_timeout, 1,
                                     JOB_TIMEOUT_LIMIT),
            failure_ttl=_bounded_int(env, 'FAILURE_TTL', cls.failure_ttl, 0, 86400),
            rate_burst=_bounded_int(env, 'RATE_BURST', cls.rate_burst, 1, 1000),
            rate_seconds=_bounded_int(env, 'RATE_SECONDS', cls.rate_seconds, 1, 3600),
            table_cache_rows=_bounded_int(env, 'TABLE_CACHE_ROWS',
                                          cls.table_cache_rows, 0, 10000000),
            table_cache_frames=_bounded_int(env, 'TABLE_CACHE_FRAMES',
                                            cls.table_cache_frames, 1, 256),
            trusted_proxy=_flag(env, 'TRUSTED_PROXY', cls.trusted_proxy),
        )
