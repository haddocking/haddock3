"""
Test general implementation in haddock3 modules.

Ensures all modules follow the same compatible architecture.
"""
import importlib.util
from pathlib import Path

import pytest

from haddock.gear.config_reader import read_config


examples_path = Path(
    Path(__file__).resolve().parents[1],
    'examples',
    )


examples_cfg_files = list(examples_path.rglob('*.cfg'))
examples_cfg_files_test = list(examples_path.rglob('*-test.cfg'))


def test_there_are_config_examples():
    """Test there are configuration files for examples."""
    assert examples_cfg_files


def test_there_are_test_config_examples():
    """Test there are configuration files for examples."""
    assert examples_cfg_files_test


@pytest.fixture(params=examples_cfg_files)
def example_config(request):
    """Give Examples's configuration files."""
    return request.param


def test_config_reader_can_read_example_configs(example_config):
    """Test gear.config_reader can read examples' configuration file."""
    read_config(example_config)


def test_integration_examples():
    """Test if integration helper scripts have the same examples."""
    run_tests = Path(examples_path, "run_tests.py")
    spec1 = importlib.util.spec_from_file_location("runtests", run_tests)
    mod1 = importlib.util.module_from_spec(spec1)
    spec1.loader.exec_module(mod1)

    compare_runs = Path(examples_path, "compare_runs.py")
    spec2 = importlib.util.spec_from_file_location("compare", compare_runs)
    mod2 = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(mod2)

    assert mod1.examples == mod2.examples
