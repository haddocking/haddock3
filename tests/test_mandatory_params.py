"""Test mandatory parameters."""
from pathlib import Path

import pytest

from haddock import haddock3_source_path
from haddock.core.defaults import RUNDIR
from haddock.gear.yaml2cfg import read_from_yaml_config
# this import is also a test
from haddock.gear.parameters import (
    MANDATORY_YAML,
    config_mandatory_general_parameters,
    mandatory_parameters,
    )


MANDATORY_YAML = Path(haddock3_source_path, "core", "mandatory.yaml")


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
