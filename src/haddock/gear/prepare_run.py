"""Logic pertraining to preparing the run files and folders."""
import logging
import shutil
from contextlib import contextmanager
from functools import wraps
from pathlib import Path

from haddock.gear.config_reader import get_module_name, read_config
from haddock.gear.parameters import config_mandatory_general_parameters
from haddock.error import ConfigurationError
from haddock.modules import modules_category
from haddock.libs.libutil import (
    copy_files_to_dir,
    make_list_if_string,
    remove_folder,
    remove_dict_keys,
    )


logger = logging.getLogger(__name__)


@contextmanager
def config_key_error():
    """Raise ConfigurationError on KeyError."""
    try:
        yield
    except KeyError as err:
        msg = f"Expected {err.args[0]!r} parameter in configuration file."
        logger.debug(err)
        raise ConfigurationError(msg) from err


def with_config_error(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        with config_key_error():
            return func(*args, **kwargs)
    return wrapper


def setup_run(workflow_path):
    """
    Setup HADDOCK3 run.

    This function performs several actions in a pipeline.

    #1 : validate the parameter TOML file
    #2 : convert strings to paths where it should
    #3 : remove folder from previous runs if run folder name overlaps
    #4 : create the needed folders/files to start the run

    Returns
    -------
    dict
        The updated parameter file.
    """
    params = read_config(workflow_path)

    validate_params(params)
    convert_params_to_path(params)
    remove_folder(params['run_dir'])
    begin_dir, _ = create_begin_files(params)

    copy_molecules_to_topology(params)

    # get a dictionary without the general config keys
    modules_params = remove_dict_keys(
        params,
        config_mandatory_general_parameters,
        )

    copy_ambig_files(modules_params, begin_dir)

    # return the modules' parameters and other parameters that may serve
    # the workflow, the "other parameters" can be expanded in the future
    # by a function if needed
    return modules_params, {'run_dir': params['run_dir']}


def validate_params(params):
    """
    Validate the parameter file.

    #1 : checks for mandatory parameters
    #2 : checks for correct modules
    """
    check_mandatory_argments_are_present(params)
    validate_modules(params)


def check_mandatory_argments_are_present(params):
    """Confirms order key exists in config."""
    for arg in config_mandatory_general_parameters:
        if arg not in params:
            _msg = (
                f"Parameter {arg!r} is not defined in the configuration file. "
                "Please refer to DOCUMENTATION-LINK for more information."
                )
            raise ConfigurationError(_msg)
    return


@with_config_error
def validate_modules(params):
    """
    Validate modules.

    Confirm the modules specified in the `order` parameter actually
    exist in HADDOCK3.

    Raises ConfigurationError if module does not exist.
    """
    keys = set(params) - set(config_mandatory_general_parameters)
    for module in keys:
        if get_module_name(module) not in modules_category.keys():
            _msg = (
                f"Module {module} not found in HADDOCK3 library. "
                "Please refer to the list of available modules at: "
                "DOCUMENTATION-LINK"
                )
            raise ConfigurationError(_msg)


def convert_params_to_path(params):
    """Convert parameters to path."""
    convert_molecules_to_path(params)
    convert_run_dir_to_path(params)
    return


@with_config_error
def convert_molecules_to_path(params):
    """
    Convert molecules path strings to Python Paths.

    And... convert `molecules` in `params` to a dictionary where keys
    are `key` + `sep` + enumerate(`start`), and values are the new Path
    values.
    """
    molecules = make_list_if_string(params['molecules'])
    params['molecules'] = [Path(i).resolve() for i in molecules]
    return


@with_config_error
def convert_run_dir_to_path(params):
    """Convert run directory value to Python Path."""
    project_dir = Path(params['run_dir'])
    params['run_dir'] = project_dir.resolve()
    return


@with_config_error
def create_begin_files(params):
    """Create initial files for HADDOCK3 run."""
    run_dir = params['run_dir']
    data_dir = run_dir / 'data'
    begin_dir = run_dir / 'begin'

    run_dir.mkdir()
    begin_dir.mkdir()
    data_dir.mkdir()

    copy_files_to_dir(params['molecules'], data_dir)
    copy_molecules_to_begin_folder(params['molecules'], begin_dir)

    return begin_dir, data_dir


def copy_molecules_to_begin_folder(
        molecules,
        begin_dir,
        mol='mol',
        sep='_',
        start=1,
        ):
    """Copy molecules to run directory begin folder."""
    for i, mol_path in enumerate(molecules, start=start):
        mol_id = f"{mol}{sep}{i}.pdb"
        begin_mol = Path(begin_dir, mol_id).resolve()
        shutil.copy(mol_path, begin_mol)


@with_config_error
def copy_molecules_to_topology(params):
    """Copy molecules to mandatory topology module."""
    params['topoaa']['molecules'] = params['molecules']


@with_config_error
def copy_ambig_files(module_params, directory):
    """Copy ambiguity table files to run directory and updates new path."""
    for step, step_dict in module_params.items():
        for key, value in step_dict.items():
            if key == 'ambig':
                ambig_f = Path(value).resolve()
                new_loc = Path(directory , step, 'ambig.tbl')
                new_loc.parent.mkdir(exist_ok=True)

                try:
                    shutil.copy(ambig_f, new_loc)
                except FileNotFoundError:
                    _msg = f'Stage: {step} ambig file {ambig_f.name} not found'
                    raise ConfigurationError(_msg)

                step_dict[key] = new_loc
