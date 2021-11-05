"""Logic pertraining to preparing the run files and folders."""
import importlib
import logging
import shutil
from contextlib import contextmanager
from functools import wraps
from pathlib import Path

from haddock import contact_us, haddock3_source_path
from haddock.core.exceptions import ConfigurationError, ModuleError
from haddock.gear.config_reader import get_module_name, read_config
from haddock.gear.parameters import config_mandatory_general_parameters
from haddock.gear.restart_run import remove_folders_after_number
from haddock.modules import (
    general_parameters_affecting_modules,
    modules_category,
    )
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


def setup_run(workflow_path, restart_from=None):
    """
    Setup HADDOCK3 run.

    This function performs several actions in a pipeline.

    #1 : validate the parameter TOML file
    #2 : convert strings to paths where it should
    #3 : copy molecules to topology key
    #4 : validate haddock3 modules params names against defaults
    #5 : remove folder from previous runs if run folder name overlaps
    #6 : create the needed folders/files to start the run
    #7 : copy additional files to run folder

    Parameters
    ----------
    workflow_path : str or pathlib.Path
        The path to the configuration file.

    erase_previous : bool
        Whether to erase the previous run folder and reprare from
        scratch. Defaults to `True`.

    Returns
    -------
    tuple of two dicts
        A dictionary with the parameters for the haddock3 modules.
        A dictionary with the general run parameters.
    """
    params = read_config(workflow_path)

    # validates the configuration file
    validate_params(params)

    # pre-treats the configuration file
    convert_params_to_path(params)
    copy_molecules_to_topology(params)

    # get a dictionary without the general config keys
    general_params = remove_dict_keys(
        params,
        list(modules_category.keys()),
        )

    modules_params = remove_dict_keys(
        params,
        list(general_params.keys()),
        )


    validate_modules_params(modules_params)
    validate_installed_modules(modules_params)

    if restart_from is None:
        # prepares the run folders
        remove_folder(general_params['run_dir'])
        begin_dir, _ = create_begin_files(general_params)

        # prepare other files
        copy_ambig_files(modules_params, begin_dir)

    else:
        remove_folders_after_number(general_params['run_dir'], restart_from)

    # return the modules' parameters and other parameters that may serve
    # the workflow, the "other parameters" can be expanded in the future
    # by a function if needed

    return modules_params, general_params


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
    keys = \
        set(params) \
        - set(config_mandatory_general_parameters) \
        - set(general_parameters_affecting_modules)

    for module in keys:
        if get_module_name(module) not in modules_category.keys():
            _msg = (
                f"Module {module} not found in HADDOCK3 library. "
                "Please refer to the list of available modules at: "
                "DOCUMENTATION-LINK"
                )
            raise ConfigurationError(_msg)


@with_config_error
def validate_modules_params(modules_params):
    """Validates individual parameters for each module."""

    for module_name, args in modules_params.items():
        _module_name = get_module_name(module_name)
        pdef = Path(
            haddock3_source_path,
            'modules',
            modules_category[_module_name],
            _module_name,
            'defaults.cfg',
            ).resolve()

        defaults = read_config(pdef)
        if not defaults:
            return

        diff = set(args.keys()) \
            - set(defaults.keys()) \
            - set(config_mandatory_general_parameters) \
            - general_parameters_affecting_modules

        if diff:
            _msg = (
                'The following parameters do not match any expected '
                f'parameters for module {module_name!r}: {diff}.'
                )
            raise ConfigurationError(_msg)


def validate_installed_modules(params):
    """Validate if third party-libraries are installed"""
    for module_name in params.keys():
        module_import_name = '.'.join([
            'haddock',
            'modules',
            modules_category[module_name],
            module_name,
            ])
        module_lib = importlib.import_module(module_import_name)
        try:
            module_lib.HaddockModule.confirm_installation()
        except Exception as err:
            _msg = (
                'A problem occurred while assessing if module '
                f'{module_name!r} is installed in your system. Have you '
                'installed the packages required to run this module? If '
                f'yes, write us at {contact_us!r} describing your system '
                'and the problems you are facing. If not, please install '
                'the required packages to use the module.'
                )
            raise ModuleError(_msg) from err


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
