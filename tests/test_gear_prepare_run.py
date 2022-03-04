"""Test prepare run module."""
from math import isnan

import pytest

from haddock.gear.prepare_run import (
    get_expandable_parameters,
    populate_mol_parameters,
    populate_topology_molecule_params,
    )
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.topology.topoaa import DEFAULT_CONFIG


DEFAULT_DICT = read_from_yaml_config(DEFAULT_CONFIG)


@pytest.mark.parametrize(
    "inp,expected",
    [
        (
            {
                "autohis": None,
                "mol1": {"nhisd", "hisd_1", "hisd_2", "nhise", "hise_1"},
                },
            {"hisd_1", "hisd_2", "hise_1"},
            )
        ]
    )
def test_get_expandable_parameters_topoaa(inp, expected):
    """Test get blocks."""
    result = get_expandable_parameters(inp, DEFAULT_DICT, "topoaa", 20)
    assert result == expected


def test_populate_topoaa_molecules():
    """Test mols are polated."""
    topoaa = {
        "molecules": ["file1.pdb", "file2.pdb"],
        "mol1": {"cyclicpept": True},
        }
    populate_topology_molecule_params(topoaa)
    assert "mol2" in topoaa
    assert topoaa["mol2"]["prot_segid"] == "B"
    assert topoaa["mol1"]["prot_segid"] == "A"
    assert topoaa["mol2"]["cyclicpept"] is False
    assert topoaa["mol1"]["cyclicpept"] is True
    assert isnan(topoaa["mol2"]["hisd_1"])
    assert isnan(topoaa["mol1"]["hisd_1"])
    assert isnan(topoaa["mol2"]["hise_1"])
    assert isnan(topoaa["mol1"]["hise_1"])
    assert topoaa["mol2"]["nhise"] == 0
    assert topoaa["mol1"]["nhise"] == 0
    assert topoaa["mol2"]["nhisd"] == 0
    assert topoaa["mol1"]["nhisd"] == 0


def test_populate_topoaa_molecules_2():
    """Test mols are polated."""
    topoaa = {
        "molecules": ["file1.pdb", "file2.pdb"],
        "mol2": {"cyclicpept": True, "prot_segid": "D"},
        }
    populate_topology_molecule_params(topoaa)
    assert "mol1" in topoaa
    assert topoaa["mol1"]["prot_segid"] == "A"
    assert topoaa["mol2"]["prot_segid"] == "D"

    assert topoaa["mol1"]["cyclicpept"] is False
    assert topoaa["mol2"]["cyclicpept"] is True

    assert isnan(topoaa["mol1"]["hisd_1"])
    assert isnan(topoaa["mol2"]["hisd_1"])
    assert isnan(topoaa["mol1"]["hise_1"])
    assert isnan(topoaa["mol2"]["hise_1"])

    assert topoaa["mol2"]["nhise"] == 0
    assert topoaa["mol1"]["nhise"] == 0
    assert topoaa["mol2"]["nhisd"] == 0
    assert topoaa["mol1"]["nhisd"] == 0


def test_populate_topoaa_molecules_3():
    """Test mols are polated."""
    topoaa = {
        "molecules": ["file1.pdb", "file2.pdb", "file3.pdb"],
        "mol2": {"cyclicpept": True, "prot_segid": "C"},
        }
    populate_topology_molecule_params(topoaa)
    assert "mol1" in topoaa
    assert topoaa["mol1"]["prot_segid"] == "A"
    assert topoaa["mol2"]["prot_segid"] == "C"
    assert topoaa["mol3"]["prot_segid"] == "B"


def test_populate_topoaa_molecules_4():
    """Test mols are polated with prot_segid sequence."""
    topoaa = {
        "molecules": ["file1.pdb", "file2.pdb", "file3.pdb", "file4.pdb"],
        "mol3": {"cyclicpept": True, "prot_segid": "A"},
        }
    populate_topology_molecule_params(topoaa)
    assert "mol1" in topoaa
    assert topoaa["mol1"]["prot_segid"] == "B"
    assert topoaa["mol2"]["prot_segid"] == "C"
    assert topoaa["mol3"]["prot_segid"] == "A"
    assert topoaa["mol4"]["prot_segid"] == "D"


def test_populate_mol_params():
    """Test populate mol."""
    params = {
        "topoaa": {"molecules": ["file1.pdb", "file2.pdb", "file3.pdb"]},
        "flexref": {"mol_fix_origin_2": True},
        "caprieval": {},
        }

    populate_mol_parameters(params)
    assert "mol_fix_origin_1" in params["flexref"]
    assert "mol_fix_origin_2" in params["flexref"]
    assert "mol_fix_origin_3" in params["flexref"]
    assert not ("mol_fix_origin_4" in params["flexref"])
    assert params["flexref"]["mol_fix_origin_2"] is True
    assert "mol_shape_1" in params["flexref"]
    assert "mol_shape_2" in params["flexref"]
    assert "mol_shape_3" in params["flexref"]
    assert not ("mol_shape_4" in params["flexref"])
    assert not params["caprieval"]
