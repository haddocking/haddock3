"""Specific tests for topoaa."""
import os
import shutil
import tempfile
from math import isnan
from pathlib import Path

import pytest

from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.topology.topoaa import DEFAULT_CONFIG as topoaa_params
from haddock.modules.topology.topoaa import HaddockModule, generate_topology

from . import golden_data


DEFAULT_DICT = read_from_yaml_config(topoaa_params)


@pytest.mark.parametrize(
    "param",
    ["hisd_1", "hise_1"],
    )
def test_variable_defaults_are_nan_in_mol1(param):
    """Test some variable defaults are as expected."""
    assert isnan(DEFAULT_DICT["mol1"][param])


def test_there_is_only_one_mol():
    """Test there is only one mol parameter in topoaa."""
    r = set(
        p
        for p in DEFAULT_DICT
        if p.startswith("mol") and p[3].isdigit()
        )
    assert len(r) == 1


@pytest.fixture
def ensemble_header_w_md5():
    return Path(golden_data, "ens_header.pdb")


@pytest.fixture
def protein():
    return Path(golden_data, "protein.pdb")


@pytest.fixture
def topoaa():
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield HaddockModule(
            order=1,
            path=Path(tmpdirname),
            initial_params=topoaa_params)


@pytest.mark.skip(reason="Not implemented yet")
def test_confirm_installation(topoaa):
    # topoaa.confirm_installation()
    pass


def test_generate_topology(topoaa, protein):
    """Test generate_topology function."""
    observed_inp_out = generate_topology(
        input_pdb=protein,
        recipe_str=topoaa.recipe_str,
        defaults=topoaa.params,
        mol_params=topoaa.params.pop('mol1'),
        default_params_path=None)

    assert observed_inp_out == Path(protein.name).with_suffix(".inp")

    observed_inp_out.unlink()


def test_get_md5(topoaa, ensemble_header_w_md5, protein):
    """Test get_md5 method."""
    observed_md5_dic = topoaa.get_md5(ensemble_header_w_md5)
    expected_md5_dic = {
        1: '71098743056e0b95fbfafff690703761',
        2: 'f7ab0b7c751adf44de0f25f53cfee50b',
        3: '41e028d8d28b8d97148dc5e548672142',
        4: '761cb5da81d83971c2aae2f0b857ca1e',
        5: '6c438f941cec7c6dc092c8e48e5b1c10'}

    assert observed_md5_dic == expected_md5_dic

    observed_md5_dic = topoaa.get_md5(protein)
    assert observed_md5_dic == {}


@pytest.mark.skip(reason="Cannot test in Github Actions")
def test__run(topoaa, protein):
    """Test _run method."""
    shutil.copy(protein, topoaa.path)
    input_mol = Path(topoaa.path, protein.name)
    current_pwd = Path.cwd()
    os.chdir(topoaa.path)
    topoaa.params['molecules'] = [input_mol]
    topoaa.envvars = topoaa.default_envvars()
    topoaa._run()
    os.chdir(current_pwd)

    inp = input_mol.with_suffix(".inp")
    out = input_mol.with_suffix(".out")
    topology = Path(input_mol.parent, f"{input_mol.stem}_haddock.psf")
    structure = Path(input_mol.parent, f"{input_mol.stem}_haddock.pdb")

    assert inp.exists()
    assert out.exists()
    assert topology.exists()
    assert structure.exists()
