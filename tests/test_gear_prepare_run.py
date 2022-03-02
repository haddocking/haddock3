"""Test prepare run module."""
import pytest

from haddock.gear.prepare_run import get_expandable_parameters
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.topology.topoaa import DEFAULT_CONFIG


@pytest.mark.parametrize(
    "inp,expected",
    [
        (
            {
                "autohis": None,
                "mol1": {"nhisd", "hisd_1", "hisd_2", "nhise", "hise_1"},
                },
            {"hisd_1", "hisd_2", "hise_1"},
            )
        ]
    )
def test_get_expandable_parameters_topoaa(inp, expected):
    """Test get blocks."""
    default_config = read_from_yaml_config(DEFAULT_CONFIG)
    result = get_expandable_parameters(inp, default_config, "topoaa")
    assert result == expected
