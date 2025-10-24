"""Test the clustfcc module."""

import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import ModuleIO
from haddock.modules.analysis.clustfcc import DEFAULT_CONFIG as clustfcc_pars
from haddock.modules.analysis.clustfcc import HaddockModule as ClustFCCModule


@pytest.fixture(name="fcc_module")
def fixture_fcc_module(monkeypatch):
    """Clustfcc module."""
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        yield ClustFCCModule(
            order=1,
            path=Path("."),
            initial_params=clustfcc_pars,
            )


def test_io_json(fcc_module, protprot_input_list):
    """Test the correct creation of the io.json file."""
    # set the input and output models
    fcc_module.previous_io.output = protprot_input_list
    fcc_module.output_models = protprot_input_list

    # export models
    fcc_module.export_io_models()
    expected_io = Path(f"{fcc_module.path}/io.json")

    assert expected_io.exists()

    # check the content of io.json
    io = ModuleIO()
    io.load(expected_io)
    assert io.input[0].file_name == protprot_input_list[0].file_name
    assert io.output[1].file_name == protprot_input_list[1].file_name


def test_one_model_one_cluster(fcc_module, protprot_input_list):
    """Test the creation of a single cluster when one input model provided."""
    single_model_list = [protprot_input_list[0]]
    # set the input and output models
    fcc_module.previous_io.output = single_model_list
    fcc_module.output_models = single_model_list
    fcc_module.run()
    assert fcc_module.output_models[0].clt_id == 1
