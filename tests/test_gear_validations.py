"""Test validation functions."""
import pytest
import random
import tempfile

from haddock.core.defaults import valid_run_dir_chars, MODULE_DEFAULT_YAML
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import Generator, Any
from haddock.gear.validations import (
    v_rundir,
    validate_parameter_scheme,
    YAML_PARAMETER_BASE_SCHEME,
    )
from haddock.libs.libio import read_from_yaml, check_yaml_duplicated_parameters


@pytest.fixture
def valid_rundir() -> str:
    """Generate a valid rundir name."""
    return ''.join(random.sample(valid_run_dir_chars, 10))


@pytest.fixture
def wrong_rundir(valid_rundir: str) -> str:
    """Generate a wrong rundir name."""
    return valid_rundir + random.choice("',!@$#%&*()|><{}[]")


@pytest.fixture
def temp_dir() -> Generator[str, Any, Any]:
    """Generate an empty directory with default.yaml."""
    with tempfile.TemporaryDirectory() as tempdir:
        yield str(tempdir)


@pytest.fixture
def dflt_duplicate_param(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with duplicated parameter."""
    fpath = f"{temp_dir}_duplicated_params_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
param1:
  key: value
param2:
  key: value
param1:
  key: value
"""
            )
    yield fpath


@pytest.fixture
def valid_default_scheme(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate valid default.yaml."""
    fpath = f"{temp_dir}_valid_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
float_param:
  default: 9.0
  min: 4.0
  max: 1000.0
  type: float
  precision: 1
  title: Floating parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
choice_string_param:
  default: 'Choice1'
  type: string
  minchars: 3
  maxchars: 10
  choices:
    - Choice1
    - Choice2
  title: Choice parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
bool_param:
  default: false
  type: boolean
  title: Boolean parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
integer_param:
  default: 10
  type: integer
  min: 1
  max: 100
  title: Integer parameter.
  short: short description.
  long: long description
  group: analysis
  explevel: expert
list_param:
  default: []
  type: list
  minitems: 0
  maxitems: 100
  title: List parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: expert
dict_param:
  dict_list_param:
    default: []
    type: list
    minitems: 0
    maxitems: 100
    title: List parameter within dict param.
    short: short descritpion
    long: long description
    group: analysis
    explevel: expert
  dict_linteger_param:
    default: 10
    type: integer
    min: 1
    max: 100
    title: Integer parameter within dict param.
    short: short description.
    long: long description
    group: analysis
    explevel: expert
  type: dict
  title: Dict parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_float_scheme1(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing float scheme max."""
    fpath = f"{temp_dir}_w_float_1_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
float_param:
  default: 9.0
  min: 4.0
  type: float
  precision: 1
  title: Floating parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
"""
            )
    yield fpath


@pytest.fixture
def w_float_scheme2(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing float scheme min."""
    fpath = f"{temp_dir}_w_float_2_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
float_param:
  default: 9.0
  max: 1000.0
  type: float
  precision: 1
  title: Floating parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
"""
            )
    yield fpath


@pytest.fixture
def w_float_scheme3(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing float scheme precision."""
    fpath = f"{temp_dir}_w_float_3_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
float_param:
  default: 9.0
  min: 4.0
  max: 1000.0
  type: float
  title: Floating parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
"""
            )
    yield fpath


@pytest.fixture
def w_int_scheme1(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing int scheme min."""
    fpath = f"{temp_dir}_w_int_1_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
integer_param:
  default: 10
  type: integer
  max: 100
  title: Integer parameter.
  short: short description.
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_int_scheme2(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing int scheme max."""
    fpath = f"{temp_dir}_w_int_2_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
integer_param:
  default: 10
  type: integer
  min: 1
  title: Integer parameter.
  short: short description.
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_base_scheme1(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme default."""
    fpath = f"{temp_dir}_w_base_1_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  type: boolean
  title: Integer parameter.
  short: short description.
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_base_scheme2(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme type."""
    fpath = f"{temp_dir}_w_base_2_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  default: false
  title: Integer parameter.
  short: short description.
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_base_scheme3(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme title."""
    fpath = f"{temp_dir}_w_base_3_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  default: false
  type: boolean
  short: short description.
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_base_scheme4(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme short."""
    fpath = f"{temp_dir}_w_base_4_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  default: false
  type: boolean
  title: Integer parameter.
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_base_scheme5(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme long."""
    fpath = f"{temp_dir}_w_base_5_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  default: false
  type: boolean
  title: Integer parameter.
  short: short description.
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_base_scheme6(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme group."""
    fpath = f"{temp_dir}_w_base_6_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  default: false
  type: boolean
  title: Integer parameter.
  short: short description.
  long: long description
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_base_scheme7(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme explevel."""
    fpath = f"{temp_dir}_w_base_7_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  default: false
  type: boolean
  title: Integer parameter.
  short: short description.
  long: long description
  group: analysis
"""
            )
    yield fpath


@pytest.fixture
def w_list_scheme1(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing list scheme minitems."""
    fpath = f"{temp_dir}_w_list_1_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
list_param:
  default: []
  type: list
  maxitems: 100
  title: List parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_list_scheme2(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing list scheme maxitems."""
    fpath = f"{temp_dir}_w_list_2_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
list_param:
  default: []
  type: list
  minitems: 0
  title: List parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: expert
"""
            )
    yield fpath


@pytest.fixture
def w_string_scheme1(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing string scheme minchars."""
    fpath = f"{temp_dir}_w_string_1_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
choice_string_param:
  default: 'Choice1'
  type: string
  maxchars: 10
  title: Choice parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
"""
            )
    yield fpath


@pytest.fixture
def w_string_scheme2(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing string scheme maxchars."""
    fpath = f"{temp_dir}_w_string_2_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
choice_string_param:
  default: 'Choice1'
  type: string
  minchars: 3
  title: Choice parameter
  short: short descritpion
  long: long description
  group: analysis
  explevel: easy
"""
            )
    yield fpath


@pytest.fixture
def w_expert_level(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with wrong expert level."""
    fpath = f"{temp_dir}_w_expert_level_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
bool_param:
  default: false
  type: boolean
  title: Integer parameter.
  short: short description.
  long: long description
  explevel: noob
"""
            )
    yield fpath


@pytest.fixture
def w_param_type(temp_dir: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with wrong parameter type."""
    fpath = f"{temp_dir}_w_param_type_{MODULE_DEFAULT_YAML}"
    with open(fpath, 'w') as filout:
        filout.write(
            """
crazy_type_param:
  default: "crazy_type"
  type: crazy_type
  title: crazy_type parameter.
  short: short description.
  long: long description
  explevel: easy
  group: crazy_param
"""
            )
    yield fpath


def test_valid_rundir(valid_rundir: str):
    """Test that valid rundir are accepted."""
    assert v_rundir(valid_rundir) is None


def test_wrong_rundir(wrong_rundir: str):
    """Test that wrong rundir throwing errors."""
    with pytest.raises(ConfigurationError):
        assert v_rundir(wrong_rundir) is None


def test_default_is_last():
    """Test that 'default' is the last check in YAML_PARAMETER_BASE_SCHEME."""
    assert YAML_PARAMETER_BASE_SCHEME[-1] == "default"


def test_check_yaml_duplicated_parameters(dflt_duplicate_param: str):
    """Test that duplicated parameters in defaults.yaml throw errors."""
    with pytest.raises(AssertionError):
        assert check_yaml_duplicated_parameters(dflt_duplicate_param) is None


def test_validate_parameter_scheme(
        valid_default_scheme: str,
        w_float_scheme1: str,
        w_float_scheme2: str,
        w_float_scheme3: str,
        w_int_scheme1: str,
        w_int_scheme2: str,
        w_string_scheme1: str,
        w_string_scheme2: str,
        w_list_scheme1: str,
        w_list_scheme2: str,
        w_base_scheme1: str,
        w_base_scheme2: str,
        w_base_scheme3: str,
        w_base_scheme4: str,
        w_base_scheme5: str,
        w_base_scheme6: str,
        w_base_scheme7: str,
        w_expert_level: str,
        w_param_type: str,
        ):
    """Test validate_parameter_scheme."""
    v_yaml_params = read_from_yaml(valid_default_scheme)
    for param_name, v_params in v_yaml_params.items():
        assert validate_parameter_scheme(param_name, v_params) is None
    all_tests_wrong_defaults = (
        w_float_scheme1,
        w_float_scheme2,
        w_float_scheme3,
        w_int_scheme1,
        w_int_scheme2,
        w_string_scheme1,
        w_string_scheme2,
        w_list_scheme1,
        w_list_scheme2,
        w_base_scheme1,
        w_base_scheme2,
        w_base_scheme3,
        w_base_scheme4,
        w_base_scheme5,
        w_base_scheme6,
        w_base_scheme7,
        w_expert_level,
        w_param_type,
        )
    for w_default_yaml_path in all_tests_wrong_defaults:
        w_yaml_params = read_from_yaml(w_default_yaml_path)
        for param_name, w_params in w_yaml_params.items():
            with pytest.raises(AssertionError):
                assert validate_parameter_scheme(param_name, w_params) is None
