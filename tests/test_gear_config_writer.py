"""Test config writer gear."""
import os
from pathlib import Path

import pytest

from haddock import EmptyPath
from haddock.gear.config_writer import (
    _convert_value_to_config_string,
    _is_dict,
    _list_by_value,
    convert_config,
    save_config,
    )


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
param7 = "string"'''


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
