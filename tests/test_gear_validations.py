"""Test validation functions."""
import pytest
import random
import tempfile

from haddock.core.defaults import valid_run_dir_chars
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import Generator, Any
from haddock.gear.validations import (
    v_rundir,
    validate_yaml_params_scheme,
    validate_defaults_yaml,
    YAML_PARAMETER_BASE_SCHEME,
    )
from haddock.libs.libio import check_yaml_duplicated_parameters


@pytest.fixture
def valid_rundir() -> str:
    """Generate a valid rundir name."""
    return ''.join(random.sample(valid_run_dir_chars, 10))


@pytest.fixture
def wrong_rundir(valid_rundir: str) -> str:
    """Generate a wrong rundir name."""
    return valid_rundir + random.choice("',!@$#%&*()|><{}[]")


@pytest.fixture
def temp_default_yaml() -> Generator[str, Any, Any]:
    """Generate an empty directory with default.yaml."""
    with tempfile.TemporaryDirectory('.') as tempdir:
        yield tempdir + "default.yaml"


@pytest.fixture
def dflt_duplicate_param(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with duplicated parameter."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def valid_default_scheme(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate valid default.yaml."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_float_scheme1(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing float scheme max."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_float_scheme2(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing float scheme min."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_float_scheme3(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing float scheme precision."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_int_scheme1(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing int scheme min."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_int_scheme2(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing int scheme max."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_base_scheme1(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme default."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_base_scheme2(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme type."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_base_scheme3(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme title."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_base_scheme4(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme short."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_base_scheme5(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme long."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_base_scheme6(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme group."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_base_scheme7(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing base scheme explevel."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_list_scheme1(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing list scheme minitems."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_list_scheme2(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing list scheme maxitems."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_string_scheme1(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing string scheme minchars."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_string_scheme2(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with missing string scheme maxchars."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_expert_level(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with wrong expert level."""
    with open(temp_default_yaml, 'w') as filout:
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
    yield temp_default_yaml


@pytest.fixture
def w_param_type(temp_default_yaml: str) -> Generator[str, Any, Any]:
    """Generate default.yaml with wrong parameter type."""
    with open(temp_default_yaml, 'w') as filout:
        filout.write(
            """
crazy_type_param:
  default: "crazy_type"
  type: crazy_type
  title: crazy_type parameter.
  short: short description.
  long: long description
  explevel: easy
"""
            )
    yield temp_default_yaml


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


def test_valid_default(valid_default_scheme: str):
    """Test that valid defaults.yaml do not throw errors."""
    assert validate_defaults_yaml(valid_default_scheme) is None


def test_wrong_float_schemes(
        w_float_scheme1: str,
        w_float_scheme2: str,
        w_float_scheme3: str,
        ):
    """Test that wrong float schemes are throwing errors."""
    with pytest.raises(ConfigurationError):
        for default in (w_float_scheme1, w_float_scheme2, w_float_scheme3):
            assert validate_defaults_yaml(default) is None
            with pytest.raises(AssertionError):
                assert validate_yaml_params_scheme(default) is None


def test_wrong_integer_schemes(
        w_int_scheme1: str,
        w_int_scheme2: str,
        ):
    """Test that wrong integer schemes are throwing errors."""
    with pytest.raises(ConfigurationError):
        for default in (w_int_scheme1, w_int_scheme2):
            assert validate_defaults_yaml(default) is None
            with pytest.raises(AssertionError):
                assert validate_yaml_params_scheme(default) is None


def test_wrong_string_schemes(
        w_string_scheme1: str,
        w_string_scheme2: str,
        ):
    """Test that wrong string schemes are throwing errors."""
    with pytest.raises(ConfigurationError):
        for default in (w_string_scheme1, w_string_scheme2):
            assert validate_defaults_yaml(default) is None
            with pytest.raises(AssertionError):
                assert validate_yaml_params_scheme(default) is None


def test_wrong_list_schemes(
        w_list_scheme1: str,
        w_list_scheme2: str,
        ):
    """Test that wrong list schemes are throwing errors."""
    with pytest.raises(ConfigurationError):
        for default in (w_list_scheme1, w_list_scheme2):
            assert validate_defaults_yaml(default) is None
            with pytest.raises(AssertionError):
                assert validate_yaml_params_scheme(default) is None


def test_wrong_base_schemes(
        w_base_scheme1: str,
        w_base_scheme2: str,
        w_base_scheme3: str,
        w_base_scheme4: str,
        w_base_scheme5: str,
        w_base_scheme6: str,
        w_base_scheme7: str,
        ):
    """Test that missing base schemes are throwing errors."""
    with pytest.raises(ConfigurationError):
        for default in (
                w_base_scheme1,
                w_base_scheme2,
                w_base_scheme3,
                w_base_scheme4,
                w_base_scheme5,
                w_base_scheme6,
                w_base_scheme7,
                ):
            assert validate_defaults_yaml(default) is None
            with pytest.raises(AssertionError):
                assert validate_yaml_params_scheme(default) is None


def test_wrong_expert_level(w_expert_level: str):
    """Test that wrong expert level is throwing errors."""
    with pytest.raises(ConfigurationError):
        assert validate_defaults_yaml(w_expert_level) is None
        with pytest.raises(AssertionError):
            assert validate_yaml_params_scheme(w_expert_level) is None


def test_param_type(w_param_type: str):
    """Test that wrong parameter type is throwing errors."""
    with pytest.raises(ConfigurationError):
        assert validate_defaults_yaml(w_param_type) is None
        with pytest.raises(AssertionError):
            assert validate_yaml_params_scheme(w_param_type) is None
