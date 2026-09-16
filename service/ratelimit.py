"""A per-client token bucket for the one route that can start work.

Design doc 3 accepts that a GET can begin a 43-second job. Until the M2 WAF exists
this is the only thing standing between an unauthenticated URL and the process pool,
so it guards the path that *starts* a job. Cache hits are cheap and stay unlimited.

The bucket table is bounded. Keying an unbounded dict on a client-controlled value
would trade a CPU exhaustion problem for a memory one.
"""
import heapq
import threading
import time

MAX_CLIENTS = 4096
# Pruning drops to this mark rather than to MAX_CLIENTS. Trimming to exactly the
# limit means the next new client prunes again, so a flood from many addresses --
# the case the limiter exists for -- would pay a full scan on every request.
LOW_WATER = MAX_CLIENTS * 3 // 4


class TokenBucket(object):
    """`burst` requests immediately, then one more every `seconds`."""

    def __init__(self, burst, seconds, clock=time.monotonic):
        self.burst = burst
        self.seconds = seconds
        self._clock = clock
        self._lock = threading.Lock()
        self._buckets = {}

    def allow(self, client):
        now = self._clock()
        with self._lock:
            tokens, last = self._buckets.get(client, (self.burst, now))
            tokens = min(self.burst, tokens + (now - last) / self.seconds)
            if tokens < 1:
                self._buckets[client] = (tokens, now)
                return False
            self._buckets[client] = (tokens - 1, now)
            if len(self._buckets) > MAX_CLIENTS:
                self._prune(now)
            return True

    def _prune(self, now):
        """Drops buckets that have refilled, so they carry no state worth keeping.

        Called with the lock held. If everything is still in use, the oldest entries
        go instead -- forgetting a bucket is fair to the client, and the alternative
        is unbounded growth.
        """
        full = self.burst * self.seconds
        for client, (_, last) in list(self._buckets.items()):
            if now - last >= full:
                del self._buckets[client]
        if len(self._buckets) > LOW_WATER:
            extra = len(self._buckets) - LOW_WATER
            oldest = heapq.nsmallest(extra, self._buckets,
                                     key=lambda c: self._buckets[c][1])
            for client in oldest:
                del self._buckets[client]


def client_ip(request, trusted_proxy=False):
    """Who to charge for this request.

    X-Forwarded-For is read only when the deployment says a proxy sets it. Behind
    the M2 ALB it carries the real client address, but with nothing in front, any
    client can set it, which would make the limiter above bypassable with a header.
    """
    if trusted_proxy:
        forwarded = request.headers.get('x-forwarded-for')
        if forwarded:
            return forwarded.split(',')[0].strip()[:64]
    return request.client.host if request.client else 'unknown'
