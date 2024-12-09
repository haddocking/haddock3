"""Test lib PDB."""
import pytest
import tempfile

from pathlib import Path

from . import data_folder
from haddock.libs.libpdb import (
    add_TER_on_chain_breaks,
    read_chainids,
    read_segids,
    )


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
    result = read_chainids(lines)
    assert result == expected


@pytest.mark.parametrize(
    "lines,expected",
    (
        (chainC, ["C", "C", "C"]),
        ),
    )
def test_read_seg_ids(lines, expected):
    result = read_segids(lines)
    assert result == expected


def test_add_TER_on_chain_breaks(monkeypatch):
    """Test proper detection of chain breaks."""
    with tempfile.TemporaryDirectory(".") as tdir:
        output_fpath = Path(tdir, "test_withTER.pdb")
        add_TER_on_chain_breaks(
            Path(data_folder, "noTER.pdb"),
            output_fpath,
            )
        # Make sure file was created
        assert output_fpath.exists()
        assert output_fpath.stat().st_size != 0
        # Make sure the file resemble the reference one
        with open(output_fpath, "r") as testfile, \
                open(Path(data_folder, "withTER.pdb")) as reffile:
            # Line by line comparisons
            for test_, ref_ in zip(testfile, reffile):
                assert test_.rstrip() == ref_.rstrip()
