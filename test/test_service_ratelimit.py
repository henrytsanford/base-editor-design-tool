"""The token bucket guarding the path that starts a job."""
from service.ratelimit import MAX_CLIENTS, TokenBucket, client_ip


class FakeRequest(object):
    def __init__(self, host='1.2.3.4', headers=None):
        self.client = type('C', (), {'host': host})()
        self.headers = headers or {}


def bucket(burst=3, seconds=10):
    clock = [0.0]
    return TokenBucket(burst, seconds, clock=lambda: clock[0]), clock


def test_a_burst_is_allowed_then_refused():
    b, _ = bucket()
    assert [b.allow('a') for _ in range(5)] == [True, True, True, False, False]


def test_tokens_come_back_over_time():
    b, clock = bucket()
    for _ in range(3):
        b.allow('a')
    assert b.allow('a') is False
    clock[0] = 10.0
    assert b.allow('a') is True


def test_clients_are_independent():
    b, _ = bucket()
    for _ in range(3):
        b.allow('a')
    assert b.allow('a') is False
    assert b.allow('b') is True


def test_the_bucket_table_stays_bounded():
    """Keyed on a client-controlled value, so it must not grow without limit."""
    b, clock = bucket(burst=1, seconds=1)
    for i in range(MAX_CLIENTS + 500):
        clock[0] += 1
        b.allow('client-%d' % i)
    assert len(b._buckets) <= MAX_CLIENTS


def test_forwarded_for_is_ignored_unless_proxies_are_declared():
    """Without a proxy in front, the header is attacker-controlled, so trusting it
    would make the limiter bypassable."""
    request = FakeRequest('1.2.3.4', {'x-forwarded-for': '9.9.9.9'})
    assert client_ip(request, proxy_hops=0) == '1.2.3.4'
    assert client_ip(request, proxy_hops=1) == '9.9.9.9'


def test_forwarded_for_is_read_from_the_right():
    """Google's front end appends to a client-supplied header, so everything left
    of the entry it wrote is whatever the caller sent."""
    spoofed = FakeRequest('10.0.0.1', {'x-forwarded-for': '6.6.6.6, 9.9.9.9'})
    assert client_ip(spoofed, proxy_hops=1) == '9.9.9.9'


def test_rotating_a_spoofed_entry_does_not_change_the_client():
    """The bypass the right-hand read exists to close: a fresh bucket per request."""
    seen = {client_ip(FakeRequest('10.0.0.1',
                                  {'x-forwarded-for': '%d.0.0.1, 9.9.9.9' % i}),
                      proxy_hops=1)
            for i in range(20)}
    assert seen == {'9.9.9.9'}


def test_two_proxies_skip_the_one_nearest_the_app():
    """A load balancer in front of Cloud Run appends its own address after the
    client's, so the client is second from the right."""
    request = FakeRequest('10.0.0.1',
                          {'x-forwarded-for': '6.6.6.6, 9.9.9.9, 130.211.0.1'})
    assert client_ip(request, proxy_hops=2) == '9.9.9.9'


def test_a_header_shorter_than_the_proxy_chain_is_not_trusted():
    """Too few entries means the expected proxies did not write it."""
    request = FakeRequest('1.2.3.4', {'x-forwarded-for': '9.9.9.9'})
    assert client_ip(request, proxy_hops=2) == '1.2.3.4'
    blank = FakeRequest('1.2.3.4', {'x-forwarded-for': '9.9.9.9, '})
    assert client_ip(blank, proxy_hops=1) == '1.2.3.4'


def test_forwarded_for_is_length_capped():
    long = FakeRequest('1.2.3.4', {'x-forwarded-for': 'a' * 500})
    assert len(client_ip(long, proxy_hops=1)) <= 64


def test_a_request_without_a_client_still_has_a_key():
    request = FakeRequest()
    request.client = None
    assert client_ip(request) == 'unknown'
