"""
Test general implementation in haddock3 modules.

Ensures all modules follow the same compatible architecture.
"""
import importlib

import pytest

from haddock import modules_defaults_path
from haddock.core.exceptions import ConfigurationError
from haddock.gear.config_reader import read_config
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules import _not_valid_config, config_readers, modules_category


@pytest.fixture(params=modules_category.items())
def module(request):
    """Give imported HADDOCK3 modules."""
    module_name, category = request.param
    mod = ".".join(['haddock', 'modules', category, module_name])
    module = importlib.import_module(mod)
    return module


def test_config_reader_can_read_defaults(module):
    """Test gear.config_reader can read modules' default file."""
    if module.DEFAULT_CONFIG.parent.name == 'topocg':
        read_from_yaml_config(module.DEFAULT_CONFIG)
    else:
        assert read_from_yaml_config(module.DEFAULT_CONFIG)


def test_general_config():
    """Test general config is readable."""
    assert read_config(modules_defaults_path)


def test_all_defaults_have_the_same_name(module):
    """Test all default configuration files have the same name."""
    assert module.DEFAULT_CONFIG.name == 'defaults.yml'


def test_config_readers_keys():
    """Test config readers have all options."""
    assert set(config_readers.keys()) == {".yml", ".cfg"}


def test_not_valid_config():
    """Test not valid config error."""
    with pytest.raises(ConfigurationError):
        with _not_valid_config():
            config_readers[".zzz"]
