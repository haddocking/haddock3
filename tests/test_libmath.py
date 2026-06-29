"""Test libmath."""

import pytest

from haddock.libs.libmath import RandomNumberGenerator


@pytest.mark.parametrize(
    "lower,upper",
    [
        (0, 9),
        (100, 99999),  # the range used for CNS seeds in libcns
        (5, 5),  # degenerate single-value range
        (-3, 3),
    ],
)
def test_randint_stays_within_bounds(lower, upper):
    """randint must never return a value outside [lower, upper]."""
    rng = RandomNumberGenerator()
    for _ in range(10000):
        value = rng.randint(lower, upper)
        assert lower <= value <= upper


def test_randint_is_reproducible_for_same_seed():
    """Same seed must yield the same sequence."""
    rng_a = RandomNumberGenerator(seed=42)
    rng_b = RandomNumberGenerator(seed=42)
    seq_a = [rng_a.randint(100, 99999) for _ in range(50)]
    seq_b = [rng_b.randint(100, 99999) for _ in range(50)]
    assert seq_a == seq_b


def test_randint_covers_full_range():
    """Both endpoints must be reachable across enough samples."""
    rng = RandomNumberGenerator()
    values = {rng.randint(0, 3) for _ in range(2000)}
    assert values == {0, 1, 2, 3}
