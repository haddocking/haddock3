"""Test yaml2cfg gear."""
import filecmp
import os
import pytest
from pathlib import Path
from fnmatch import fnmatch

from haddock import haddock3_source_path
from haddock.libs.libio import read_from_yaml
from haddock.gear.yaml2cfg import (
    flat_yaml_cfg,
    yaml2cfg_text,
    )
from haddock.gear.validations import validate_defaults_yaml

from . import (
    haddock3_yaml_cfg_examples,
    haddock3_yaml_converted,
    haddock3_yaml_converted_no_header,
    )


@pytest.fixture
def default_yaml_files():
    """Return list of defaults.yaml file within the haddock src directory."""
    all_defaults_yaml: list[str] = []
    default_yaml_fname = "defaults.yaml"
    for path, _subdirs, files in os.walk(haddock3_source_path):
        for name in files:
            if fnmatch(name, default_yaml_fname):
                all_defaults_yaml.append(f"{path}/{default_yaml_fname}")
    return all_defaults_yaml


complex_cfg = {
    "param1": {
        "default": 1,
        "other": None,
        "others": [None, None],
        "explevel": "easy",
        },
    "param2": {
        "default": 2,
        "other": None,
        "others": [None, None],
        "explevel": "easy",
        },
    "param3": {
        "param4": {
            "default": 4,
            "other": None,
            "others": [None, None],
            "explevel": "easy",
            },
        },
    "param5": {
        "param6": {
            "param7": {
                "default": 7,
                "other": None,
                "others": [None, None],
                "explevel": "easy",
                },
            },
        },
    }


complex_cfg_simplified = {
    "param1": 1,
    "param2": 2,
    "param3": {"param4": 4},
    "param5": {"param6": {"param7": 7}},
    }


no_explvl_key = {
    "param1": {
        "default": 1,
        "other": None,
        "others": [None, None],
        },
    }


def test_flat_complex_config_1():
    """Test if complex config is flatten properly."""
    result = flat_yaml_cfg(complex_cfg)
    assert result == complex_cfg_simplified


def test_yaml2cfg_test():
    """Test yaml dict to cfg."""
    ycfg = read_from_yaml(haddock3_yaml_cfg_examples)
    result = yaml2cfg_text(ycfg, "topoaa", "all")
    assert isinstance(result, str)

    p = Path('dummy_test.cfg')
    p.write_text(result)

    assert filecmp.cmp(p, haddock3_yaml_converted, shallow=False)
    p.unlink()


def test_yaml2cfg_test_no_header():
    """Test yaml dict to cfg."""
    ycfg = read_from_yaml(haddock3_yaml_cfg_examples)
    result = yaml2cfg_text(ycfg, None, "all")
    assert isinstance(result, str)

    p = Path('dummy_test_2.cfg')
    p.write_text(result)

    assert filecmp.cmp(
        p,
        haddock3_yaml_converted_no_header,
        shallow=False,
        )
    p.unlink()


def test_yaml_duplicated_params(default_yaml_files):
    """Make sure no duplicated parameters are present in a ymal file."""
    assert default_yaml_files != []
    # Loop over default yaml files
    for default_yaml_fpath in default_yaml_files:
        validate_defaults_yaml(default_yaml_fpath)
