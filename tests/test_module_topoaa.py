"""Specific tests for topoaa."""
from math import isnan

import pytest

from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.topology.topoaa import DEFAULT_CONFIG


DEFAULT_DICT = read_from_yaml_config(DEFAULT_CONFIG)


@pytest.mark.parametrize(
    "param",
    ["hisd_1", "hise_1"],
    )
def test_variable_defaults_are_nan_in_mol1(param):
    """Test some variable defaults are as expected."""
    assert isnan(DEFAULT_DICT["mol1"][param])


def test_there_is_only_one_mol():
    """Test there is only one mol parameter in topoaa."""
    r = set(p for p in DEFAULT_DICT if p.startswith("mol"))
    assert len(r) == 1
