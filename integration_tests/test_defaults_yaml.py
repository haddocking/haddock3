"""Test contents of all defaults.yaml module parameters."""
import os
import pytest
import tempfile
from fnmatch import fnmatch

from haddock import haddock3_source_path
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import Generator, Any
from haddock.gear.validations import (
    validate_defaults_yaml,
    validate_yaml_params_scheme,
    )


@pytest.fixture
def default_yaml_files() -> list[str]:
    """Return list of defaults.yaml file within the haddock src directory."""
    all_defaults_yaml: list[str] = []
    # Loop over all the files in src/haddock/
    for path, _, files in os.walk(haddock3_source_path):
        for name in files:
            # Check if it is a defaults.yaml file
            if fnmatch(name, MODULE_DEFAULT_YAML):
                all_defaults_yaml.append(f"{path}/{MODULE_DEFAULT_YAML}")
    return all_defaults_yaml


@pytest.fixture
def temp_dir() -> Generator[str, Any, Any]:
    """Generate an empty directory with default.yaml."""
    with tempfile.TemporaryDirectory() as tempdir:
        yield str(tempdir)


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


def test_yaml_duplicated_params(default_yaml_files: list[str]):
    """Make sure no duplicated parameters are present in a ymal file."""
    assert default_yaml_files != []
    # Loop over default yaml files
    for default_yaml_fpath in default_yaml_files:
        assert validate_defaults_yaml(default_yaml_fpath) is None


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
        for default in (w_float_scheme1, w_float_scheme2, w_float_scheme3, ):
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
