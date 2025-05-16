"""Specific tests for topocg."""

import os
import tempfile
from math import isnan
from pathlib import Path
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import pytest

from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.topology.topocg import DEFAULT_CONFIG as topocg_params
from haddock.modules.topology.topocg import HaddockModule as Topocg
from haddock.modules.topology.topocg import generate_topology

from haddock.libs.libontology import Format

from . import golden_data


DEFAULT_DICT = read_from_yaml_config(topocg_params)


@pytest.mark.parametrize(
    "param",
    ["hisd_1", "hise_1"],
)
def test_variable_defaults_are_nan_in_mol1(param):
    """Test some variable defaults are as expected."""
    assert isnan(DEFAULT_DICT["mol1"][param])


def test_there_is_only_one_mol():
    """Test there is only one mol parameter in topocg."""
    r = set(p for p in DEFAULT_DICT if p.startswith("mol") and p[3].isdigit())
    assert len(r) == 1


@pytest.fixture(name="ensemble_header_w_md5")
def fixture_ensemble_header_w_md5():
    """???"""
    return Path(golden_data, "ens_header.pdb")


@pytest.fixture(name="protein")
def fixture_protein():
    """???"""
    return Path(golden_data, "protein.pdb")


@pytest.fixture(name="topocg")
def fixture_topocg(monkeypatch):
    """topocg module fixture"""
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        yield Topocg(order=1, path=Path("."), initial_params=topocg_params)


def test_generate_topology(topocg, protein):
    """Test generate_topology function."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        force_field=topocg.params["cgffversion"]
        observed_inp_out = generate_topology(
            input_pdb=protein,
            output_path=topocg.path,
            recipe_str=topocg.recipe_str,
            defaults=topocg.params,
            mol_params=topocg.params.pop("mol1"),
            default_params_path=None,
            force_field=force_field,
        )

    assert observed_inp_out == Path(protein.name).with_suffix(".inp")
    
    # Confirm the expected output file exists
    expected_filename = topocg.path.resolve() /".." / f"{protein.stem}_cg.pdb" #Path(f"/tmp/{protein.stem}_cg.pdb")
    assert expected_filename.exists(), f"Expected CG PDB file not found: {expected_filename}"

    # Inspect the content to confirm it's a CG structure
    with open(expected_filename, "r") as f:
        contents = f.read()

    # Check for CG-specific atom/residue names (e.g., "BB" for backbone bead)
    assert "BB" in contents, "CG markers not found in output PDB."


def test_get_md5(topocg, ensemble_header_w_md5, protein):
    """Test get_md5 method."""
    observed_md5_dic = topocg.get_md5(ensemble_header_w_md5)
    expected_md5_dic = {
        1: "71098743056e0b95fbfafff690703761",
        2: "f7ab0b7c751adf44de0f25f53cfee50b",
        3: "41e028d8d28b8d97148dc5e548672142",
        4: "761cb5da81d83971c2aae2f0b857ca1e",
        5: "6c438f941cec7c6dc092c8e48e5b1c10",
    }

    assert observed_md5_dic == expected_md5_dic

    observed_md5_dic = topocg.get_md5(protein)
    assert observed_md5_dic == {}
