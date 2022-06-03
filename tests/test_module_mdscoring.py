"""Test the mdscoring module."""
import os
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.mdscoring import DEFAULT_CONFIG as mdscoring_pars
from haddock.modules.scoring.mdscoring import HaddockModule

from . import golden_data


@pytest.fixture
def output_models():
    """Prot-prot models using for mdscoring output."""
    return [
        PDBFile(
            Path(golden_data, "protprot_complex_1.pdb"),
            path=golden_data,
            score=42.0
            ),
        PDBFile(
            Path(golden_data, "protprot_complex_2.pdb"),
            path=golden_data,
            score=28.0
            )]


def test_mdscoring_output(output_models):
    """Test emscoring expected output."""
    mds_module = HaddockModule(
        order=1,
        path=Path("1_emscoring"),
        initial_params=mdscoring_pars
        )
    # original names
    mds_module.output_models = output_models
    for mod in range(len(output_models)):
        ori_name = "original_name_" + str(mod) + ".pdb"
        mds_module.output_models[mod].ori_name = ori_name
    # creating output
    output_fname = Path("mdscoring.tsv")
    mds_module.output(output_fname)
    observed_outf_l = [e.split() for e in open(
        output_fname).readlines() if not e.startswith('#')]
    # expected output
    expected_outf_l = [
        ["structure", "original_name", "md5", "score"],
        ["protprot_complex_1.pdb", "original_name_0.pdb", "None", "42.0"],
        ["protprot_complex_2.pdb", "original_name_1.pdb", "None", "28.0"]]

    assert observed_outf_l == expected_outf_l

    os.unlink(output_fname)
