"""Test HaddockModel gear."""

from pathlib import Path

import pytest

from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libontology import PDBFile

from . import golden_data


@pytest.fixture
def protprot_input_list():
    """Prot-prot input."""
    return [
        PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data),
    ]


@pytest.fixture
def protprot_1bkd_input_list():
    """
    Prot-prot input for target 1bkd.

    Heterogeneous ensemble and big protein.
    """
    return [
        PDBFile(Path(golden_data, "protprot_1bkd_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protprot_1bkd_2.pdb"), path=golden_data),
    ]


def test_haddockmodel_energies(protprot_input_list):
    """Test haddock model."""
    haddock_mod = HaddockModel(protprot_input_list[0].rel_path)
    assert isinstance(haddock_mod.energies, dict)

    exp_energies = {
        "vdw": 1.85709,
        "elec": -9.41558,
        "desolv": 3.25569,
        "air": 482.494,
        "coup": 0.0,
        "sym": 0.0,
        "rdcs": 0.0,
        "vean": 0.0,
        "dani": 0.0,
        "xpcs": 0.0,
        "rg": 0.0,
        "bsa": 846.821,
        "total": 474.936,
        "bonds": 0.0,
        "angles": 0.0,
        "improper": 0.0,
        "dihe": 0.0,
        "cdih": 0.0,
    }

    assert haddock_mod.energies == exp_energies


def test_nonhaddockmodel(protprot_1bkd_input_list):
    """Test a non haddock model."""
    haddock_mod = HaddockModel(protprot_1bkd_input_list[0].rel_path)
    assert haddock_mod.energies == {}


def test_haddockmodel_score(protprot_input_list):
    """Test haddock model score."""
    haddock_mod = HaddockModel(protprot_input_list[1].rel_path)

    weights = {
        "w_vdw": 1.0,
        "w_elec": 1.0,
        "w_desolv": 1.0,
        "w_air": 1.0,
    }

    assert haddock_mod.calc_haddock_score(**weights) == 84.9855

    weights["w_air"] = 0.1
    weights["w_bsa"] = -0.01

    assert haddock_mod.calc_haddock_score(**weights) == -13.38146
