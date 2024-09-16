"""Test the clustfcc module."""

import os
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import ModuleIO, PDBFile
from haddock.modules.analysis.clustfcc import DEFAULT_CONFIG as clustfcc_pars
from haddock.modules.analysis.clustfcc import HaddockModule as ClustFCCModule

from . import golden_data


@pytest.fixture
def prot_input_list():
    """Prot-prot input."""
    return [
        PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data),
    ]


@pytest.fixture
def fcc_module():
    """Clustfcc module."""
    with tempfile.TemporaryDirectory() as tempdir:
        yield ClustFCCModule(order=1, path=Path(tempdir), initial_params=clustfcc_pars)


def test_io_json(fcc_module, prot_input_list):
    """Test the correct creation of the io.json file."""
    # set the input and output models
    fcc_module.previous_io.output = prot_input_list
    fcc_module.output_models = prot_input_list

    # export models
    fcc_module.export_io_models()

    assert Path("io.json").exists()

    # check the content of io.json
    io = ModuleIO()
    io.load("io.json")
    assert io.input[0].file_name == prot_input_list[0].file_name
    assert io.output[1].file_name == prot_input_list[1].file_name

    # remove io.json file
    os.unlink(Path("io.json"))
