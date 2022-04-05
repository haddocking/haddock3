"""Test libio."""
import pytest

from haddock.libs.libio import read_from_yaml

from . import emptycfg, haddock3_yaml_cfg_examples


@pytest.mark.parametrize(
    "cfg",
    [
        emptycfg,
        haddock3_yaml_cfg_examples,
        ],
    )
def test_read_from_yaml(cfg):
    """Test read from yaml file."""
    result = read_from_yaml(cfg)
    assert isinstance(result, dict)
