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


def test_forwarded_for_is_ignored_unless_a_proxy_is_trusted():
    """Without a proxy in front, the header is attacker-controlled, so trusting it
    would make the limiter bypassable."""
    request = FakeRequest('1.2.3.4', {'x-forwarded-for': '9.9.9.9'})
    assert client_ip(request, trusted_proxy=False) == '1.2.3.4'
    assert client_ip(request, trusted_proxy=True) == '9.9.9.9'


def test_forwarded_for_takes_the_first_hop_and_is_length_capped():
    request = FakeRequest('1.2.3.4', {'x-forwarded-for': '9.9.9.9, 10.0.0.1'})
    assert client_ip(request, trusted_proxy=True) == '9.9.9.9'
    long = FakeRequest('1.2.3.4', {'x-forwarded-for': 'a' * 500})
    assert len(client_ip(long, trusted_proxy=True)) <= 64


def test_a_request_without_a_client_still_has_a_key():
    request = FakeRequest()
    request.client = None
    assert client_ip(request) == 'unknown'
