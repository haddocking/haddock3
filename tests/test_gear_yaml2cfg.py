"""Test yaml2cfg gear."""
import filecmp
import re
import os
import pytest
from pathlib import Path
from fnmatch import fnmatch

from haddock.gear.yaml2cfg import flat_yaml_cfg, yaml2cfg_text
from haddock.libs.libio import read_from_yaml
from haddock import haddock3_source_path

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
    # Build regular expression
    yaml_param_regex = re.compile("^(([A-Za-z0-9]_?)+):")
    for yaml_fpath in default_yaml_files:
        # Loop over default yaml files
        parsed_param_names: dict[str, int] = {}
        with open(yaml_fpath, 'r') as filin:
            yaml_content = filin.readlines()
        for i, line in enumerate(yaml_content, start=1):
            if (match := yaml_param_regex.search(line)):
                param_name = match.group(1)
                assert param_name not in parsed_param_names.keys(), f"Parameter '{param_name}' in {yaml_fpath} has duplicates: l.{parsed_param_names[param_name]} and l.{i}"  # noqa : E501
                parsed_param_names[param_name] = i
