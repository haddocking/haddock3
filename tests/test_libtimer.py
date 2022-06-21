"""Test libtimer."""
import pytest

from haddock.libs.libtimer import convert_seconds_to_min_sec, log_time


def test_logtime():
    with log_time("some string"):
        2 * 2


@pytest.mark.parametrize(
    "seconds,expected",
    [
        (60, "1 minute and 0 seconds"),
        (120, "2 minutes and 0 seconds"),
        (40, "40 seconds"),
        (179, "2 minutes and 59 seconds"),
        (3600, "1h0m0s"),
        (3601, "1h0m1s"),
        (3600 + 120, "1h2m0s"),
        (3600 + 125, "1h2m5s"),
        (3600 * 2 + 125, "2h2m5s"),
        (3600 * 2 + 179, "2h2m59s"),
        (3600 * 2 + 180, "2h3m0s"),
        ]
    )
def test_convert_seconds(seconds, expected):
    """Convert seconds to min&sec."""
    result = convert_seconds_to_min_sec(seconds)
    assert result == expected
