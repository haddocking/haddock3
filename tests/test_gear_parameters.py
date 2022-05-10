"""Test mandatory parameters."""
import pytest

from haddock.core.defaults import RUNDIR
# this import is also a test
from haddock.gear.parameters import (  # noqa: F401
    MANDATORY_YAML,
    OPTIONAL_YAML,
    _mandatory_parameters,
    _optional_parameters,
    config_mandatory_general_parameters,
    config_optional_general_parameters,
    )
from haddock.gear.yaml2cfg import read_from_yaml_config


def test_mandatory_config_set():
    """Test if is set."""
    assert isinstance(config_mandatory_general_parameters, set)


def test_optional_config_set():
    """Test if is set."""
    assert isinstance(config_optional_general_parameters, set)


def test_rundir_in_yaml():
    """Test RUNDIR is in yaml."""
    assert RUNDIR in _mandatory_parameters


def test_rundir_in_set():
    """Test RUNDIR is in yaml."""
    assert RUNDIR in config_mandatory_general_parameters


def test_optional_skip_preprocessing_yaml():
    assert "skip_preprocess" in _optional_parameters


def test_optional_skip_preprocessing_set():
    assert "skip_preprocess" in config_optional_general_parameters
