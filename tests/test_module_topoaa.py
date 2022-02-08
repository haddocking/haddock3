"""Specific tests for topoaa."""
from math import isnan

import pytest

from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.topology.topoaa import DEFAULT_CONFIG


DEFAULT_DICT = read_from_yaml_config(DEFAULT_CONFIG)


@pytest.mark.parametrize(
    "idx",
    list(map(str, range(1, 21))),
    )
def test_variable_defaults_are_nan(idx):
    """Test some variable defaults are as expected."""
    assert isnan(DEFAULT_DICT[f"mol{idx}"]["hisd_1"])
    assert isnan(DEFAULT_DICT[f"mol{idx}"]["hise_1"])
