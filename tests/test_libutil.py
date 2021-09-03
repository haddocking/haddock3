import pytest

from haddock.libs.libutil import non_negative_int


@pytest.mark.parametrize(
    'i,expected',
    [
        (1, 1),
        ('1', 1),
        (0, 0),
        ('0', 0),
        ],
    )
def test_non_negative_int(i, expected):
    """."""
    r = non_negative_int(i)
    assert r == expected


@pytest.mark.parametrize('i', [-1, -12, -14])
def test_non_negative_int_error(i):
    """."""
    with pytest.raises(ValueError):
        non_negative_int(i)
