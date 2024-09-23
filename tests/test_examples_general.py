"""
Test general implementation in haddock3 modules.

Ensures all modules follow the same compatible architecture.
"""
import importlib.util
import pytest

from pathlib import Path

from haddock.gear.config import load as read_config
from haddock.gear.prepare_run import (
    ALL_POSSIBLE_GENERAL_PARAMETERS,
    identify_modules,
    validate_module_names_are_not_misspelled,
    validate_parameters_are_not_misspelled,
    validate_modules_params,
    )
from haddock.gear.parameters import config_optional_general_parameters_dict
from haddock.libs.libutil import recursive_dict_update, remove_dict_keys
from haddock.modules import non_mandatory_general_parameters_defaults


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


@pytest.fixture(params=examples_cfg_files, name="example_config")
def fixture_example_config(request):
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
    assert mod1.examples == mod2.examples, f"{mod1.examples}!={mod2.examples}"


def test_validate_cfg_files():
    """Test all the examples configuration files are valid."""
    for cfg_file in examples_cfg_files:
        assert validate_cfg_file(cfg_file), f"Error detected in {cfg_file}!"


def validate_cfg_file(cfg_file: Path) -> bool:
    """Handeler to validate a configuration file.

    Parameters
    ----------
    cfg_file : Path
        Path to the configuration file to check.

    Returns
    -------
    bool
        True if config file is OK else False
    """
    try:
        if cfg_file.name == "params.cfg":
            # skip the params.cfg files as they might be part of old runs
            return True
        config_files = read_config(cfg_file)
        # update default non-mandatory parameters with user params
        params = recursive_dict_update(
            config_optional_general_parameters_dict,
            config_files["final_cfg"],
            )

        params = recursive_dict_update(
            non_mandatory_general_parameters_defaults,
            params,
            )

        validate_module_names_are_not_misspelled(params)

        # separate general from modules' parameters
        _modules_keys = identify_modules(params)
        general_params = remove_dict_keys(params, _modules_keys)
        modules_params = remove_dict_keys(params, list(general_params.keys()))

        validate_parameters_are_not_misspelled(
            general_params,
            reference_parameters=ALL_POSSIBLE_GENERAL_PARAMETERS,
            )

        validate_modules_params(modules_params, 20)
    except Exception as e:
        print(e)
        return False
    else:
        return True
