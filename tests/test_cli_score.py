"""Tests related to haddock.clis.cli_score"""
import pytest

from haddock.clis.cli_score import get_parameters
from haddock.core.exceptions import ConfigurationError
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.scoring.emscoring import (
    DEFAULT_CONFIG as EMSCORING_DEFAULTS_CONFIG_PATH,
    )


@pytest.fixture
def empty_params() -> dict:
    return {}


@pytest.fixture
def v_cmd_line_params() -> dict:
    return {
        "w_bsa": 10,
        "w_desolv": 10,
        "w_elec": 10,
        "w_vdw": 10,
        }


@pytest.fixture
def wrong_params() -> dict:
    return {"fake": "wrong"}


@pytest.fixture
def wrong_params_type() -> dict:
    return {"w_bsa": "wrong"}


@pytest.fixture
def default_emscoring_params() -> dict:
    default_emscoring = read_from_yaml_config(EMSCORING_DEFAULTS_CONFIG_PATH)
    return default_emscoring


def test_no_input_params(empty_params: dict, default_emscoring_params: dict):
    """Test get_parameters without inputs."""
    final_params = get_parameters(empty_params)
    assert isinstance(final_params, dict)
    for param_name, param_value in default_emscoring_params.items():
        assert final_params[param_name] == param_value


def test_input_params(
        v_cmd_line_params: dict[str, int],
        default_emscoring_params: dict,
        ):
    """Test get_parameters with inputs."""
    final_params = get_parameters(v_cmd_line_params)
    assert isinstance(final_params, dict)
    for param_name, param_value in default_emscoring_params.items():
        if param_name in v_cmd_line_params.keys():
            assert final_params[param_name] == v_cmd_line_params[param_name]
        else:
            assert final_params[param_name] == param_value


def test_wrong_params(wrong_params: dict[str, str]):
    """Test get_parameters with wrong inputs."""
    with pytest.raises(ConfigurationError):
        final_params = get_parameters(wrong_params)
        assert final_params is None


def test_wrong_param_type(wrong_params_type: dict[str, str]):
    """Test get_parameters with wrong inputs."""
    with pytest.raises(ConfigurationError):
        final_params = get_parameters(wrong_params_type)
        assert final_params is None
