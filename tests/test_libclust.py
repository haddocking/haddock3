"""Test the libclust library."""
import os
from pathlib import Path

import pytest

from haddock.libs.libclust import write_unclustered_list
from haddock.libs.libontology import PDBFile

from . import golden_data


@pytest.fixture
def protprot_input_models():
    """Prot-prot input."""
    return [
        PDBFile(
            Path(golden_data, "protprot_complex_1.pdb"),
            path=golden_data,
            score=-42),
        PDBFile(
            Path(golden_data, "protprot_complex_2.pdb"),
            path=golden_data,
            score=-17
            )
        ]


@pytest.fixture
def protprot_output_models():
    """Prot-prot output."""
    return [
        PDBFile(
            Path(golden_data, "protprot_complex_1.pdb"),
            path=golden_data,
            score=-42
            )
        ]


def test_write_unclustered_list(protprot_input_models, protprot_output_models):
    """Test write_unclustered_list function."""
    write_unclustered_list(protprot_input_models, protprot_output_models)
    uncl_fname = "unclustered.txt"
    observed_file_content = open(uncl_fname, "r").read()
    expected_file_content = (
        f"### Unclustered structures ###{os.linesep}"
        f"{os.linesep}"
        f"1\tprotprot_complex_2.pdb\t-17.00{os.linesep}"
        f"{os.linesep}"
        )
    assert observed_file_content == expected_file_content
    os.unlink(uncl_fname)
