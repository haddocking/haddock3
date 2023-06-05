"""Test the clustfcc module."""
import os
from pathlib import Path

import pytest

# from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.clustfcc import DEFAULT_CONFIG as clustfcc_pars
from haddock.modules.analysis.clustfcc import HaddockModule

# from . import golden_data


@pytest.fixture
def protprot_input_list():
    """Prot-prot input."""
    return []


@pytest.fixture
def output_list():
    """Clustfcc output list."""
    return [
        "fcc.matrix",
        "clustfcc.txt",
        "io.json",
        "clustfcc.tsv",
    ]


@pytest.fixture
def fcc_module(protprot_input_list, output_list):
    fcc_module = HaddockModule(
        order=1, path=Path("1_emscoring"), initial_params=clustfcc_pars
    )
    fcc_module.previous_io.output = protprot_input_list

    yield fcc_module

    for f in output_list:
        path_f = Path(f)
        if path_f.exists():
            os.unlink(path_f)


@pytest.mark.skip(reason="Not implemented yet")
def test_clustfcc_module(fcc_module, output_list):
    """Test the clustfcc module."""

    fcc_module._run()

    # Check if the output is correct and all files are there
    ls = os.listdir()
    for expected_file in output_list:
        assert expected_file in ls, "Expected file not found"

    # Check the contact files are correct
    expected_output_length = [100, 119]
    observed_output_lengths = []

    con1_len = len(open("protprot_complex_1.con").readlines())
    observed_output_lengths.append(con1_len)

    con2_len = len(open("protprot_complex_2.con").readlines())
    observed_output_lengths.append(con2_len)

    assert observed_output_lengths == expected_output_length
