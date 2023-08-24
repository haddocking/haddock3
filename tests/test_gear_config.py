"""Test config reader."""
from math import isnan
from pathlib import Path

import pytest

from haddock.gear import config


@pytest.mark.parametrize(
    'line,expected',
    [
        ('[header]', 'header'),
        ('[under_scores_work]', 'under_scores_work'),
        ],
    )
def test_main_header_re(line, expected):
    """Test header regex."""
    result = config._main_header_re.match(line)
    assert result[1] == expected


@pytest.mark.parametrize(
    'line,expected,expected2',
    [
        ('[another.header]', 'another', '.header'),
        ('[another.header_2]', 'another', '.header_2'),
        ('[another.header.h3]', 'another', '.header.h3'),
        ],
    )
def test_sub_header_re(line, expected, expected2):
    """Test header regex."""
    group = config._sub_header_re.match(line)
    assert group[1] == expected
    assert group[2] == expected2


@pytest.mark.parametrize(
    'line,expected',
    [
        ("['header.3']", 'header'),
        ("['header.1234']", 'header'),
        ],
    )
def test_main_quoted_header_re(line, expected):
    """Test quoted header regex."""
    result = config._main_quoted_header_re.match(line)
    assert result[1] == expected


@pytest.mark.parametrize(
    'line,expected1,expected2',
    [
        ("['another.3'.header]", 'another', '.header'),
        ("['another.3'.header.some]", 'another', '.header.some'),
        ],
    )
def test_sub_quoted_header_re(line, expected1, expected2):
    """Test sub quoted header regex."""
    result = config._sub_quoted_header_re.match(line)
    assert result[1] == expected1
    assert result[2] == expected2


@pytest.mark.parametrize(
    'line',
    [
        '[[header]]',
        '[héàder]',
        'value = "some_string"',
        '[header with spaces]',
        '[not.valid]',
        ],
    )
def test_main_header_re_wrong(line):
    """Test main header wrong."""
    assert config._main_header_re.match(line) is None


@pytest.mark.parametrize(
    'line',
    [
        '[[header]]',
        'value = "some_string"',
        '[header with spaces]',
        '[not.valid with spaces]',
        '[single]',
        ],
    )
def test_sub_header_re_wrong(line):
    """Test sub header wrong."""
    assert config._sub_header_re.match(line) is None


@pytest.mark.parametrize(
    'header,name',
    [
        ('header', 'header'),
        ('header.1', 'header'),
        ('header.1.1..1', 'header'),
        ],
    )
def test_get_module_name(header, name):
    """Test get module name."""
    assert config.get_module_name(header) == name


_config_example_1 = """
# some comment
num1 = 10
name="somestring"
w_vdw_1 = 1.0
w_vdw_0 = 0.01
w_vdw_2 = 1.0
[headerone]
name = "theotherstring"
_list = [
    12,
    15,
    ]

[headerone.weights]
val = 1
bsa = 0.1
elec = -1
list1 = [1, 2,3]

[headerone] #some comment
[headerone.weights]
other = 50
"""

_config_example_dict_1 = {
    "num1": 10,
    "name": "somestring",
    "w_vdw_0": 0.01,
    "w_vdw_1": 1.0,
    "w_vdw_2": 1.0,
    "headerone.1": {
        "name": "theotherstring",
        "_list": [12, 15],
        "weights": {
            "val": 1,
            "bsa": 0.1,
            "elec": -1,
            "list1": [1, 2, 3],
            },
        },
    "headerone.2": {'weights': {'other': 50}},
    }

_config_example_2 = """num1 = 10
molecules = ["../some/file", "../some/otherfile"]
[module]
some_path_fname = "./pointing/to/some/path"
[module.d1]
var1 = 1
[module.d1.d2]
var3 = true
list_ = [
    1,
    2,
    3, 4
    ]
"""


_config_example_dict_2 = {
    "num1": 10,
    "molecules": [Path("../some/file"), Path("../some/otherfile")],
    "module.1": {
        "some_path_fname": Path("pointing", "to", "some", "path"),
        "d1": {
            "var1": 1,
            "d2": {
                "var3": True,
                "list_": [1, 2, 3, 4],
                },
            },
        },
    }

# this examples shows the behaviour of subkey repetition
_config_example_3 = """
val = 1
[header]
[header.d1]
val2 = 10
val3 = 20
"""

_config_example_dict_3 = {
    "val": 1,
    "header.1": {
        "d1": {
            "val2": 10,
            "val3": 20,
            }
        }
    }


@pytest.mark.parametrize(
    'config_example,expected',
    [
        (_config_example_1, _config_example_dict_1),
        (_config_example_2, _config_example_dict_2),
        (_config_example_3, _config_example_dict_3),
        ],
    )
def test_load(config_example, expected):
    """Test read config."""
    r = config.loads(config_example)
    assert r['final_cfg'] == expected


def test_load_nan_vlaue():
    """Test read config."""
    r = config.loads("param=nan")
    assert isnan(r['final_cfg']["param"])


@pytest.mark.parametrize(
    "conf",
    [
        _config_example_dict_1,
        _config_example_dict_2,
        _config_example_dict_3,
        ]
    )
def test_save_config(conf):
    target = Path("dummy_config.cfg")
    config.save(conf, target)
    loaded = config.load(target)
    assert loaded['final_cfg'] == conf
    target.unlink()


@pytest.mark.parametrize(
    "conf",
    [
        _config_example_dict_1,
        _config_example_dict_2,
        _config_example_dict_3,
        ]
    )
def test_save_config_toml(conf):
    target = Path("dummy_config.toml")
    config.save(conf, target, pure_toml=True)
    target.unlink()


@pytest.mark.parametrize(
    "paramline",
    [
        'param=true',
        'sampling_factor=false',
        'sym3 =true',
        'sym4 = false',
        ]
    )
def test_not_uppercase_bool(paramline):
    """Test lower case regex."""
    assert config._uppercase_bool_re.match(paramline) is None


@pytest.mark.parametrize(
    "paramline,expectedparam,expectedbool",
    [
        ('param=True', 'param=', 'true'),
        ('param = False', 'param = ', 'false'),
        ('sampling_factor =False', 'sampling_factor =', 'false'),
        ('sym3=True', 'sym3=', 'true'),
        ('_strange_1_param_= False', '_strange_1_param_= ', 'false'),
        ('_strange1_2_3param_ = True', '_strange1_2_3param_ = ', 'true'),
        ]
    )
def test_uppercase_bool(paramline, expectedparam, expectedbool):
    """Test lower case regex."""
    groups = config._uppercase_bool_re.match(paramline)
    assert groups is not None  # Able to catch line
    paramgroup = groups[1]
    assert paramgroup == expectedparam  # Able to catch parameter
    lowercase_bool = groups[4].lower()
    assert lowercase_bool == expectedbool  # Able to lowercase boolean
