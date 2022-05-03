"""Test config writer gear."""
import collections.abc
import importlib
import os
from math import isnan
from pathlib import Path

import pytest
import toml

from haddock import EmptyPath
# from haddock.gear.config_reader import read_config
from haddock.gear.config_writer import (
    _convert_value_to_config_string,
    _is_dict,
    _list_by_value,
    convert_config,
    save_config,
    )
from haddock.gear.yaml2cfg import read_from_yaml_config

from . import working_modules


@pytest.mark.parametrize(
    "value, expected",
    [
        (range(10), False),
        (10, False),
        (10.0, False),
        ("a", False),
        (["a", 1], False),
        (set(range(10)), False),
        (tuple(range(10)), False),
        (Path.cwd(), False),
        ({"a": None}, True),
        ]
    )
def test_is_dict(value, expected):
    result = _is_dict((None, value))
    assert result == expected


@pytest.mark.parametrize(
    "value, expected",
    [
        (False, "false"),
        (1, "1"),
        (1.9, "1.9"),
        ("param", '"param"'),
        (Path('path'), f'"{str(Path("path").resolve())}"'),
        (EmptyPath(), '""'),
        (float('nan'), "nan"),
        (None, "none"),
        ]
    )
def test_convert_value_to_config_string(value, expected):
    result = _convert_value_to_config_string(value, fmt="{}")
    assert result == expected


@pytest.mark.parametrize(
    "value",
    (
        {None: None},
        list(range(10)),
        tuple(range(10)),
        )
    )
def test_convert_value_to_config_string_AssertionError(value):
    with pytest.raises(AssertionError):
        _convert_value_to_config_string(value, fmt="{}")


case_1 = {
    "param1": 1,
    "param2": 1.0,
    "param3": False,
    "param4": ["a", "b", "c"],
    "sub1": {"param5": 50},
    "param6": float('nan'),
    "param7": "string",
    "param8": [],
    "param9": tuple([]),
    }

case_1_text = '''param1 = 1
param2 = 1.0
param3 = false
param4 = [
    "a",
    "b",
    "c"
    ]

param6 = nan
param7 = "string"
param8 = []
param9 = []

[sub1]
param5 = 50'''


def test_convert_config_1():
    """Test basic parameter config conversion."""
    result = os.linesep.join(convert_config(case_1))
    assert result == case_1_text


case_1_text_ignore = '''param1 = 1
param3 = false
param4 = [
    "a",
    "b",
    "c"
    ]

param6 = nan
param7 = "string"
param8 = []
param9 = []'''


def test_convert_config_1_ignore():
    """Test basic parameter config conversion."""
    _ = convert_config(case_1, ignore_params=["param2", "sub1"])
    result = os.linesep.join(_)
    assert result == case_1_text_ignore


case_2 = {
    "module": {
        "param1": 1,
        "param2": 1.0,
        "param3": False,
        "param4": ["a", "b", "c"],
        "sub1": {"param5": 50},
        "param6": float('nan'),
        "param7": "string",
        "param8": [],
        "param9": tuple([]),
        },
    }

case_2_text = '''[module]
param1 = 1
param2 = 1.0
param3 = false
param4 = [
    "a",
    "b",
    "c"
    ]

param6 = nan
param7 = "string"
param8 = []
param9 = []
[module.sub1]
param5 = 50'''


def test_convert_config_2():
    """Test with module."""
    result = os.linesep.join(convert_config(case_2, module_names=["module"]))
    assert result == case_2_text


values_1 = ["a", 10, Path.cwd()]
values_str = ['    "a",', "    10,", "    " + f'"{str(Path.cwd().resolve())}"']


def test_list_by_value_1():
    result = list(_list_by_value(values_1))
    assert result == values_str


def test_list_by_value_AssertionError():
    with pytest.raises(AssertionError):
        list(_list_by_value([]))


@pytest.mark.parametrize(
    "d,text, kws",
    [
        (case_1, case_1_text, {}),
        (case_2, case_2_text, {"module_names": ["module"]}),
        ]
    )
def test_save_config(d, text, kws):
    fpath = Path("save_config_test.cfg")
    save_config(d, fpath, **kws)
    result = fpath.read_text()
    assert result == text
    fpath.unlink()


def test_save_config_module_name():
    fpath = Path("save_config_test.cfg")
    save_config(case_1, fpath, module_name="module")
    result = fpath.read_text()
    assert result == case_2_text
    fpath.unlink()


def test_save_config_TypeError():
    with pytest.raises(TypeError):
        save_config(case_1, "file.cfg", module_name=10)


@pytest.fixture(params=working_modules)
def modules(request):
    """Give imported HADDOCK3 modules."""
    module_name, category = request.param
    mod = ".".join(['haddock', 'modules', category, module_name])
    module = importlib.import_module(mod)
    return module


@pytest.fixture()
def module_yaml_dict(modules):
    """Give flatted yaml config dictionaries."""
    cfg = read_from_yaml_config(modules.DEFAULT_CONFIG)
    return cfg


def test_load_save_configs(module_yaml_dict):
    """
    Test reading, writing, and reloading a config.

    Test compatibility between read_yaml_config on the defaults
    and saving config.

    Step 1: reads the default configs from the modules using
        read_from_yaml_config function. This gives a dictionary with
        built-in python types.

    Step 2: save the step 1 config dictionary to a haddock3 user config
        file using the save_config function

    Step 3: reads the config file from step 2 using `toml` library.

    Step 4: compares dictionaries from step 1 and step 3 recursively

    Notes: Why do we use `toml` to read the config file and not our
        gear.config_reader.read_config function? Because the latter converts
        strings that point to paths to Path and EmptyPath objects.

        TODO: maybe the read_from_yaml_config function should convert
        those also? Like not because of the interactions of haddock with
        third-party software. Or yes with options.
        @joaomcteixeira Apr 21, 2022.

        Hence, we use tolm so we can directly compare python types.
    """
    pcfg = Path('module_test.cfg')
    save_config(module_yaml_dict, pcfg)
    # rcfg = read_config(pcfg)
    rcfg = toml.load(pcfg)
    _assert_dict(rcfg, module_yaml_dict)
    pcfg.unlink()


def _assert_dict(d1, d2):
    """Compare dicts recursively."""
    d1keys = set(d1.keys())
    d2keys = set(d2.keys())

    assert not d1keys.difference(d2keys), "keys in configs differ"

    for key, value in d1.items():
        if isinstance(value, collections.abc.Mapping):
            _assert_dict(value, d2[key])
            continue

        _assert_dict_with_nan(d1[key], d2[key])


def _assert_dict_with_nan(d1v, d2v):
    # nan cannot be compared with "=="
    assert d1v == d2v or (isnan(d1v) and isnan(d2v))
