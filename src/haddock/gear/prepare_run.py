"""Logic pertraining to preparing the run files and folders."""
import logging
import shutil
from contextlib import contextmanager
from functools import wraps
from pathlib import Path

import toml

from haddock import modules_folder
from haddock.error import ConfigurationError
from haddock.gear import mandatory_params
from haddock.libs.libutil import copy_files_to_dir, remove_folder


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
    params = toml.load(workflow_path)

    validate_params(params)
    convert_params_to_path(params)
    remove_folder(params['run_dir'])
    create_begin_files(params)

    return params


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
    for arg in mandatory_params.general:
        if arg not in params:
            _msg = (
                f"Parameter {arg!r} is not defined in the configuration file. "
                "Please refer to DOCUMENTATION-LINK for more information."
                )
            raise ConfigurationError(_msg)
    return


@with_config_error
def validate_modules(params):
    """Validate modules."""

    modes = (
        (package, params['stage'][package].get('mode', 'default'))
        for package in params['order']
        )

    for package, mode in modes:
        mode = mode or 'default'
        module_loc = Path(modules_folder, package, mode).with_suffix('.py')
        if not module_loc.exists():
            _msg = (
                f"Method {package}:{mode} not found in HADDOCK3 library. "
                "Please refer to the list of available modules at: "
                "DOCUMENTATION-LINK"
                )
            raise ConfigurationError(_msg)
        else:
            params['stage'][package]['mode'] = mode


def convert_params_to_path(params):
    """Convert parameters to path."""
    convert_molecules_to_path(params)
    convert_run_dir_to_path(params)
    return


@with_config_error
def convert_molecules_to_path(params, key='mol', sep='_', start=1):
    """
    Convert molecules path strings to Python Paths.

    And... convert `molecules` in `params` to a dictionary where keys
    are `key` + `sep` + enumerate(`start`), and values are the new Path
    values.
    """
    new_paths = (
        (f'{key}{sep}{i}', Path(file_name))
        for i, file_name in enumerate(params['molecules'], start)
        )

    params['molecules'] = dict(new_paths)
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


    copy_files_to_dir(params['molecules'].values(), data_dir)
    copy_molecules_to_begin_folder(params['molecules'], begin_dir)
    copy_ambig_files(params, begin_dir)

    return


def copy_molecules_to_begin_folder(mol_dict, begin_dir):
    """Copy molecules to run directory."""
    for mol_id, mol_path in mol_dict.items():

        begin_mol = Path(begin_dir, f'{mol_id}.pdb').resolve()
        shutil.copy(mol_path, begin_mol)
        mol_dict[mol_id] = begin_mol


@with_config_error
def copy_ambig_files(params, directory):
    """Copy ambiguity table files to run directory and updates new path."""
    for step, step_dict in params['stage'].items():
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
