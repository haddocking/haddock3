"""Test lib PDB."""
import pytest

from haddock.libs import libpdb


chainC = [
    'ATOM      3  CA  ARG C   4      37.080  43.455  -3.421  1.00  0.00      C    C  ',  # noqa: E501
    'ATOM      3  CA  GLU C   6      33.861  45.127  -2.233  1.00  0.00      C    C  ',  # noqa: E501
    'ATOM      3  CA  ALA C   7      35.081  45.036   1.305  1.00  0.00      C    C  ',  # noqa: E501
    ]


@pytest.mark.parametrize(
    "lines,expected",
    (
        (chainC, ["C", "C", "C"]),
        ),
    )
def test_read_chain_ids(lines, expected):
    result = libpdb.read_chainids(lines)
    assert result == expected


@pytest.mark.parametrize(
    "lines,expected",
    (
        (chainC, ["C", "C", "C"]),
        ),
    )
def test_read_seg_ids(lines, expected):
    result = libpdb.read_segids(lines)
    assert result == expected
