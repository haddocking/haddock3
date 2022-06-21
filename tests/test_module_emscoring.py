"""Test the emscoring module."""
import os
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.emscoring import DEFAULT_CONFIG as emscoring_pars
from haddock.modules.scoring.emscoring import HaddockModule

from . import golden_data


@pytest.fixture
def output_models():
    """Prot-DNA models using for emscoring output."""
    return [
        PDBFile(
            Path(golden_data, "protdna_complex_1.pdb"),
            path=golden_data,
            score=42.0
            ),
        PDBFile(
            Path(golden_data, "protdna_complex_2.pdb"),
            path=golden_data,
            score=28.0
            )]


def test_emscoring_output(output_models):
    """Test emscoring expected output."""
    ems_module = HaddockModule(
        order=1,
        path=Path("1_emscoring"),
        initial_params=emscoring_pars
        )
    # original names
    ems_module.output_models = output_models
    for mod in range(len(output_models)):
        ori_name = "original_name_" + str(mod) + ".pdb"
        ems_module.output_models[mod].ori_name = ori_name
    # creating output
    output_fname = Path("emscoring.tsv")
    ems_module.output(output_fname)
    observed_outf_l = [e.split() for e in open(
        output_fname).readlines() if not e.startswith('#')]
    # expected output
    expected_outf_l = [
        ["structure", "original_name", "md5", "score"],
        ["protdna_complex_1.pdb", "original_name_0.pdb", "None", "42.0"],
        ["protdna_complex_2.pdb", "original_name_1.pdb", "None", "28.0"]]

    assert observed_outf_l == expected_outf_l

    os.unlink(output_fname)
