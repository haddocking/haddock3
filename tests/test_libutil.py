from pathlib import Path

import pytest

from haddock.libs.libutil import file_exists, non_negative_int


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


@pytest.mark.parametrize(
    'i,expected',
    [
        (Path(__file__), Path(__file__)),
        (str(Path(__file__)), Path(__file__)),
        ],
    )
def test_file_exists(i, expected):
    """."""
    r = file_exists(i)
    assert r == expected


@pytest.mark.parametrize(
    'i',
    [
        'some_bad.path',
        Path(__file__).parent, # this is a folder
        ],
    )
def test_file_exists_wrong(i):
    """."""
    with pytest.raises(ValueError):
        file_exists(i)
