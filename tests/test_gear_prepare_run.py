"""Test prepare run module."""
import shutil
from math import isnan
from pathlib import Path

import pytest

from haddock.core.exceptions import ConfigurationError
from haddock.gear.prepare_run import (
    check_if_path_exists,
    copy_molecules_to_topology,
    fuzzy_match,
    get_expandable_parameters,
    populate_mol_parameters,
    populate_topology_molecule_params,
    update_step_contents_to_step_names,
    validate_module_names_are_not_misspelled,
    validate_param_range,
    validate_param_type,
    validate_parameters_are_not_misspelled,
    validate_value,
    )
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules import modules_names
from haddock.modules.topology.topoaa import DEFAULT_CONFIG

from . import data_folder, steptmp


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
        "topoaa.1": {"molecules": ["file1.pdb", "file2.pdb", "file3.pdb"]},
        "flexref.1": {"mol_fix_origin_2": True},
        "caprieval.1": {},
        }

    populate_mol_parameters(params)
    assert "mol_fix_origin_1" in params["flexref.1"]
    assert "mol_fix_origin_2" in params["flexref.1"]
    assert "mol_fix_origin_3" in params["flexref.1"]
    assert not ("mol_fix_origin_4" in params["flexref.1"])
    assert params["flexref.1"]["mol_fix_origin_2"] is True
    assert "mol_shape_1" in params["flexref.1"]
    assert "mol_shape_2" in params["flexref.1"]
    assert "mol_shape_3" in params["flexref.1"]
    assert not ("mol_shape_4" in params["flexref.1"])
    assert not params["caprieval.1"]


def test_check_if_path_exists():
    Path("file_01.txt").write_text("a")
    Path("file_02.txt").write_text("a")
    check_if_path_exists("file_01.txt")
    check_if_path_exists("file_02.txt")

    with pytest.raises(ValueError) as err:
        check_if_path_exists("file-01.txt")
        end = ("the following \'file_01.txt\' is the closest match to the "
               "supplied \'file-01.txt\', did you mean to open this?")
        assert str(err).endswith(end)

    Path("file_01.txt").unlink()
    Path("file_02.txt").unlink()


@pytest.mark.parametrize(
    "user_input,expected",
    [
        ("long-fromat", [("long-fromat", "long-format")]),
        (["loong-format", "verboese"],
            [("loong-format", "long-format"), ("verboese", "verbose")]),
        ("middle-format", [("middle-format", "long-format")]),
        ("out", [("out", "output-dir")]),
        ]
    )
def test_fuzzy_match(user_input, expected):
    possibilities = ["long-format", "short-format", "verbose", "output-dir"]
    assert fuzzy_match(user_input, possibilities) == expected


@pytest.mark.parametrize(
    "molecules,expected",
    [
        ("mol1.pdb", 1),
        (["mol1.pdb", "mol2.pdb"], 2),
        ],
    )
def test_copy_mols_to_topo_dir(molecules, expected):
    d = {}
    copy_molecules_to_topology(molecules, d)
    assert "molecules" in d
    assert len(d["molecules"]) == expected
    assert all(isinstance(m, Path) for m in d["molecules"])


def test_validate_step_names_are_not_misspelled():
    params = {
        "par1": None,
        }

    for i, module in enumerate(modules_names):
        # the key and value for this function do not matter
        if i % 2 == 0:
            params[module] = {None: None}
        else:
            # add index
            params[module + ".1"] = {None: None}

    validate_module_names_are_not_misspelled(params)


def test_validate_step_names_are_not_misspelled_error():
    params = {
        "par1": None,
        }

    params["r1g1dbod1"] = {None: None}

    with pytest.raises(ValueError):
        validate_module_names_are_not_misspelled(params)


def test_validate_params():
    params = ["param1", "param2"]
    ref = ["param1", "param2", "param3"]
    validate_parameters_are_not_misspelled(params, ref)


def test_validate_params_error():
    params = ["param1", "param4"]
    ref = ["param1", "param2", "param3"]
    with pytest.raises(ValueError):
        validate_parameters_are_not_misspelled(params, ref)


def test_update_step_folders_from_restart():
    output_tmp = Path(data_folder, "1_dummystep")
    shutil.copytree(steptmp, output_tmp)
    prev_names = ['0_dummystep']
    next_names = ['1_dummystep']
    update_step_contents_to_step_names(prev_names, next_names, data_folder)

    file1 = Path(output_tmp, "file1.in").read_text()
    assert '0_dummystep' not in file1
    assert '1_dummystep' in file1

    file2 = Path(output_tmp, "folder", "file2.in").read_text()
    assert '0_dummystep' not in file2
    assert '1_dummystep' in file2

    shutil.rmtree(output_tmp)


@pytest.mark.parametrize(
    "defaultparam,inputvalue",
    [
        ({'type': 'integer'}, 1),
        ({'type': 'float'}, 0.99),
        ({'type': 'list'}, [1, 2]),
        ({'type': 'boolean'}, False),
        ({'type': 'boolean'}, True),
        ({'type': 'string'}, 'Hello world'),
        ]
    )
def test_accepted_value_type(defaultparam, inputvalue):
    """Test if the provided value type is correct.

    Should not return anything
    """
    assert validate_param_type(defaultparam, inputvalue) is None


@pytest.mark.parametrize(
    "defaultparam,inputvalue",
    [
        ({'type': 'integer', 'default': 0}, False),
        ({'type': 'float', 'default': 1.1}, 'Hello world'),
        ({'type': 'list', 'default': [1, 2]}, {1: 1, 2: 2}),
        ({'type': 'boolean', 'default': True}, 1),
        ({'type': 'boolean', 'default': False}, 0),
        ({'type': 'str', 'default': 'test'}, 2.1),
        ]
    )
def test_value_type_error(defaultparam, inputvalue):
    """Test if the provided value type is wrong.

    Should return a formated string describing the error.
    """
    assert validate_param_type(defaultparam, inputvalue) is not None
    

@pytest.mark.parametrize(
    "defaultparam,inputvalue",
    [
        ({'min': 0, 'max': 1}, 0.5),
        ({'min': 1, 'max': 1000}, 1),
        ({'min': 0, 'max': 1000}, 1000),
        ({'min': 0, 'max': 1000}, -0),
        ({'choices': ['ok', 'fine']}, 'fine'),
        ({'choices': ['super']}, 'super'),
        ]
    )
def test_accepted_value_range(defaultparam, inputvalue):
    """Test if the provided value range/choices is correct.

    Should not return anything
    """
    assert validate_param_range(defaultparam, inputvalue) is None


@pytest.mark.parametrize(
    "defaultparam,inputvalue",
    [
        ({'min': 0, 'max': 1, 'default': 0}, 2),
        ({'min': 1.1, 'max': 1.3, 'default': 1.2}, 1.4),
        ({'min': 0.1, 'max': 1.0, 'default': 0.5}, -0.1),
        ({'choices': ['not', 'good'], 'default': 'good'}, ['not', 'good']),
        ({'choices': ['super'], 'default': 'super'}, 'wrong'),
        ]
    )
def test_value_range_error(defaultparam, inputvalue):
    """Test if the provided value range/choices is erronnated.

    Should return formated string describing the error.
    """
    assert validate_param_range(defaultparam, inputvalue) is not None


@pytest.mark.parametrize(
    "defaultparams,key,value",
    [
        ({'key': {'min': 0, 'max': 1, 'type': 'integer'}}, 'key', 1),
        ({'key': {'min': 0.5, 'max': 2.2, 'type': 'float'}}, 'key', 1.1),
        ({'key': {'min': 0, 'max': 1, 'type': 'integer'}}, 'key', 0),
        ({'key': {'choices': ['ok', 'fine'], 'type': 'string'}},
         'key', 'ok'),
        ({'key': {'choices': ['ok', 'fine'], 'type': 'string'}},
         'key', 'fine'),
        ({'key': {'type': 'boolean'}}, 'key', True),
        ({'key': {'type': 'list'}}, 'key', ['mol1', 'mol2']),
        ({'mol': {'group': 'molecules', 'type': 'list'}}, 'mol', ['pdb']),
        ({'molecule': {'type': 'list'}}, 'molecules', ['mol1', 'pdb']),
        ]
    )
def test_accepted_param_value(defaultparams, key, value):
    """Test if the provided parameter = value is correct.

    Should not return anything
    """
    assert validate_value(defaultparams, key, value) is None


@pytest.mark.parametrize(
    "defaultparams,key,value",
    [
        ({'key': {'min': 0, 'max': 1, 'type': 'integer',
                  'default': 0}}, 'key', 1.1),
        ({'key': {'min': 0.5, 'max': 2.2,
                  'type': 'float', 'default': 1.1}}, 'key', 0),
        ({'key': {'choices': ['not', 'good'],
                  'type': 'string',
                  'default': 'good'}},
         'key', 'ok'),
        ({'key': {'choices': ['False', 'True'],
                  'type': 'string',
                  'default': 'True'}},
         'key', True),
        ({'key': {'choices': ['False', 'True'], 'type': 'string',
                  'default': 'False'}},
         'key', False),
        ({'key': {'type': 'boolean', 'default': False}}, 'key', 0),
        ({'key': {'type': 'boolean', 'default': True}}, 'key', 1),
        ({'key': {'type': 'boolean', 'default': False}}, 'key', []),
        ({'key': {'type': 'boolean', 'default': False}}, 'key', ''),
        ({'key': {'type': 'list', 'default': ['m1', 'm2']}}, 'key', 1),
        ({'key': {'type': 'list', 'default': ['m1', 'm2']}}, 'key', 'list'),
        ({'key': {'type': 'list', 'default': ['m1', 'm2']}},
         'key', {'hello': 'world'}),
        ]
    )
def test_param_value_error(defaultparams, key, value):
    """Test if the provided parameter = value is not correct.

    Should raise ConfigurationError
    """
    with pytest.raises(ConfigurationError):
        validate_value(defaultparams, key, value)
