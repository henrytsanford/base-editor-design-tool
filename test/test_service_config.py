"""Settings read from the environment."""
import pytest

from service.config import MAX_PROXY_HOPS, Settings


def test_proxy_hops_default_to_ignoring_the_header():
    assert Settings.from_env({}).trusted_proxy_hops == 0


def test_proxy_hops_are_read_and_bounded():
    assert Settings.from_env({'TRUSTED_PROXY_HOPS': '1'}).trusted_proxy_hops == 1
    with pytest.raises(ValueError):
        Settings.from_env({'TRUSTED_PROXY_HOPS': str(MAX_PROXY_HOPS + 1)})
    with pytest.raises(ValueError):
        Settings.from_env({'TRUSTED_PROXY_HOPS': '-1'})


@pytest.mark.parametrize('value', ['true', '1', 'yes', 'on', 'TRUE'])
def test_asking_to_trust_the_header_the_old_way_fails_at_startup(value):
    """A deployment that expects the header trusted must not silently run with a
    limiter that ignores it."""
    with pytest.raises(ValueError, match='TRUSTED_PROXY_HOPS'):
        Settings.from_env({'TRUSTED_PROXY': value})


@pytest.mark.parametrize('value', ['false', '0', ''])
def test_trusted_proxy_false_still_starts(value):
    """It asks for what the default already does."""
    assert Settings.from_env({'TRUSTED_PROXY': value}).trusted_proxy_hops == 0


def test_result_limits_are_read_and_bounded():
    settings = Settings.from_env({'TABLE_MAX_ROWS': '1000', 'RESULTS_MAX_MB': '64'})
    assert settings.table_max_rows == 1000 and settings.results_max_mb == 64
    with pytest.raises(ValueError):
        Settings.from_env({'RESULTS_MAX_MB': '0'})
