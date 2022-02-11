"""Test mandatory parameters."""
import pytest

from haddock.core.defaults import RUNDIR
# this import is also a test
from haddock.gear.parameters import (  # noqa: F401
    MANDATORY_YAML,
    _mandatory_parameters,
    config_mandatory_general_parameters,
    )
from haddock.gear.yaml2cfg import read_from_yaml_config


@pytest.fixture
def mandatory_parameters_dict():
    """Test mandatory params config."""
    return read_from_yaml_config(MANDATORY_YAML)


def test_config_set():
    """Test if is set."""
    assert isinstance(config_mandatory_general_parameters, set)


def test_rundir_in_yaml(mandatory_parameters_dict):
    """Test RUNDIR is in yaml."""
    assert RUNDIR in mandatory_parameters_dict


def test_rundir_in_set():
    """Test RUNDIR is in yaml."""
    assert RUNDIR in config_mandatory_general_parameters
