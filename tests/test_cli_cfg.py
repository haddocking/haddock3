"""Test haddock3-cfg client."""
import pytest

from haddock import config_expert_levels
from haddock.clis import cli_cfg
from haddock.modules import modules_category


@pytest.fixture(params=config_expert_levels + ("all",))
def config_level(request):
    """Haddock3 config levels."""
    return request.param


@pytest.fixture(params=(True, False))
def global_params(request):
    """Haddock3 config levels."""
    return request.param


@pytest.fixture(params=list(modules_category.keys()) + [None])
def module(request):
    return request.param


def test_export_cfgs_add_global(module, config_level, global_params):
    """Test export all configs work with `add_global` parameter."""
    cli_cfg.main(module, explevel=config_level, global_params=global_params)


def test_export_cfgs_details(module, config_level, global_params):
    """Test export all configs work with detailed parameter description."""
    cli_cfg.main(
        module,
        explevel=config_level,
        global_params=global_params,
        details=True,
        )
