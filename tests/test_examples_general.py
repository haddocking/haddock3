"""
Test general implementation in haddock3 modules.

Ensures all modules follow the same compatible architecture.
"""
import importlib
from pathlib import Path

import pytest

from haddock.gear.config_reader import read_config


examples_path = Path(
    Path(__file__).resolve().parents[1],
    'examples',
    )


examples_cfg_files = list(examples_path.rglob('*.cfg'))


def test_there_are_config_examples():
    """Test there are configuration files for examples."""
    assert examples_cfg_files


@pytest.fixture(params=examples_cfg_files)
def example_config(request):
    return request.param


def test_config_reader_can_read_example_configs(example_config):
    """Test gear.config_reader can read examples' configuration file."""
    read_config(example_config)
