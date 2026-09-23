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
MAX_PROXY_HOPS = 5


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
    # Results with more guides than this are offered as downloads only, never parsed
    # into a table. The default preset makes ~20x the guides of an NGG editor, and
    # TTN at edit=all is 517k of them: a 687 MB frame. 50,000 rows is ~65 MB, and
    # leaves every result of an ordinary gene viewable.
    table_max_rows: int = 50000
    # Size budget for RESULTS_DIR, which is in memory on Cloud Run. Past it, the
    # least recently viewed results are deleted after each job (storage.evict).
    results_max_mb: int = 1024
    # How many proxies in front append to X-Forwarded-For. 0, the default, ignores
    # the header: with nothing in front, any client can set it. Cloud Run alone is 1.
    # See ratelimit.client_ip.
    trusted_proxy_hops: int = 0

    @classmethod
    def from_env(cls, env=None):
        # Defaults are read off the fields above rather than repeated here, so
        # Settings() and Settings.from_env() cannot drift apart.
        env = os.environ if env is None else env
        if (env.get('TRUSTED_PROXY') or '').strip().lower() in ('1', 'true', 'yes', 'on'):
            # Refused rather than ignored. A deployment setting it expects the header
            # to be trusted, and would otherwise run with a limiter that silently
            # ignores it; better to find out at startup. 'false' asks for what
            # TRUSTED_PROXY_HOPS=0 does anyway, so it passes.
            raise ValueError('TRUSTED_PROXY is no longer read; set TRUSTED_PROXY_HOPS '
                             'to the number of proxies in front (Cloud Run alone is 1)')
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
            table_max_rows=_bounded_int(env, 'TABLE_MAX_ROWS', cls.table_max_rows,
                                        0, 10000000),
            results_max_mb=_bounded_int(env, 'RESULTS_MAX_MB', cls.results_max_mb,
                                        1, 1000000),
            trusted_proxy_hops=_bounded_int(env, 'TRUSTED_PROXY_HOPS',
                                            cls.trusted_proxy_hops, 0, MAX_PROXY_HOPS),
        )
