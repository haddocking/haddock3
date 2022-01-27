"""Logic pertraining to preparing the run files and folders."""
import importlib
import shutil
import sys
from contextlib import contextmanager
from copy import deepcopy
from functools import wraps
from pathlib import Path

from haddock import contact_us, haddock3_source_path, log
from haddock.core.defaults import RUNDIR
from haddock.core.exceptions import ConfigurationError, ModuleError
from haddock.gear.config_reader import get_module_name, read_config
from haddock.gear.greetings import get_goodbye_help
from haddock.gear.parameters import config_mandatory_general_parameters
from haddock.gear.restart_run import remove_folders_after_number
from haddock.gear.validations import v_rundir
from haddock.libs.libutil import (
    recursive_dict_update,
    remove_dict_keys,
    zero_fill,
    )
from haddock.modules import (
    modules_category,
    non_mandatory_general_parameters_defaults,
    )


@contextmanager
def config_key_error():
    """Raise ConfigurationError on KeyError."""
    try:
        yield
    except KeyError as err:
        msg = f"Expected {err.args[0]!r} parameter in configuration file."
        raise ConfigurationError(msg) from err


def with_config_error(func):
    """Add config error context."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        with config_key_error():
            return func(*args, **kwargs)
    return wrapper


def setup_run(workflow_path, restart_from=None):
    """
    Set up HADDOCK3 run.

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
    # read config
    params = read_config(workflow_path)

    check_mandatory_argments_are_present(params)
    validate_module_names_are_not_mispelled(params)
    check_specific_validations(params)

    # update default non-mandatory parameters with user params
    params = recursive_dict_update(
        non_mandatory_general_parameters_defaults,
        params)

    clean_rundir_according_to_restart(params[RUNDIR], restart_from)

    # copy molecules parameter to topology module
    copy_molecules_to_topology(params)

    # separate general from modules parameters
    _modules_keys = identify_modules(params)
    general_params = remove_dict_keys(params, _modules_keys)
    modules_params = remove_dict_keys(params, list(general_params.keys()))

    # validations
    validate_modules_params(modules_params)
    check_if_modules_are_installed(modules_params)

    # create datadir
    data_dir = create_data_dir(general_params[RUNDIR])
    new_mp = copy_input_files_to_data_dir(data_dir, modules_params)

    # return the modules' parameters and general parameters separately
    return new_mp, general_params


def validate_params(params):
    """
    Validate the parameter file.

    #1 : checks for mandatory parameters
    #2 : checks for correct modules
    """
    check_mandatory_argments_are_present(params)
    validate_modules_names(params)


def check_mandatory_argments_are_present(params):
    """Confirm order key exists in config."""
    for arg in config_mandatory_general_parameters:
        if arg not in params:
            _msg = (
                f"Parameter {arg!r} is not defined in the configuration file. "
                "Please refer to DOCUMENTATION-LINK for more information."
                )
            raise ConfigurationError(_msg)
    return


@with_config_error
def validate_modules_names(params):
    """Validate all modules names are spelled correctly."""
    keys = \
        set(params) \
        - set(config_mandatory_general_parameters) \
        - set(non_mandatory_general_parameters_defaults)

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
    """
    Validate individual parameters for each module.

    Raises
    ------
    ConfigError
        If there is any parameter given by the user that is not defined
        in the defaults.cfg of the module.
    """
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

        blocks = get_blocks(defaults)
        block_params = read_blocks(blocks, args)

        diff = set(args.keys()) \
            - set(defaults.keys()) \
            - set(config_mandatory_general_parameters) \
            - set(non_mandatory_general_parameters_defaults.keys()) \
            - block_params

        if diff:
            _msg = (
                'The following parameters do not match any expected '
                f'parameters for module {module_name!r}: {", ".join(diff)}.'
                )
            raise ConfigurationError(_msg)


def check_if_modules_are_installed(params):
    """Validate if third party-libraries are installed."""
    for module_name in params.keys():
        module_import_name = '.'.join([
            'haddock',
            'modules',
            modules_category[get_module_name(module_name)],
            get_module_name(module_name),
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


# depecrated
# def convert_params_to_path(params):
#     """Convert parameters to path."""
#     convert_molecules_to_path(params)
#     convert_run_dir_to_path(params)
#     return
#
#
# @with_config_error
# def convert_molecules_to_path(params):
#     """
#     Convert molecules path strings to Python Paths.
#
#     And... convert `molecules` in `params` to a dictionary where keys
#     are `key` + `sep` + enumerate(`start`), and values are the new Path
#     values.
#     """
#     molecules = make_list_if_string(params['molecules'])
#     params['molecules'] = [Path(i).resolve() for i in molecules]
#     return
#
#
# @with_config_error
# def convert_run_dir_to_path(params):
#     """Convert run directory value to Python Path."""
#     params[RUNDIR] = Path(params[RUNDIR])
#     return


@with_config_error
def create_data_dir(run_dir):
    """
    Create initial files for HADDOCK3 run.

    Returns
    -------
    pathlib.Path
        A path referring only to 'data'.
    """
    data_dir = Path(run_dir, 'data')
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


@with_config_error
def copy_molecules_to_topology(params):
    """Copy molecules to mandatory topology module."""
    params['topoaa']['molecules'] = list(map(Path, params['molecules']))


def copy_input_files_to_data_dir(data_dir, modules_params):
    """Copy files to data directory."""
    new_mp = deepcopy(modules_params)
    # this line must be synchronized with create_data_dir()
    rel_data_dir = data_dir.name

    for i, molecule in enumerate(modules_params['topoaa']['molecules']):
        end_path = Path(data_dir, '00_topoaa')
        end_path.mkdir(parents=True, exist_ok=True)
        name = Path(molecule).name
        shutil.copy(molecule, Path(end_path, name))
        new_mp['topoaa']['molecules'][i] = Path(rel_data_dir, '00_topoaa', name)

    # topology always starts with 0
    for i, (module, params) in enumerate(modules_params.items(), start=0):
        end_path = Path(f'{zero_fill(i)}_{get_module_name(module)}')
        for parameter, value in params.items():
            if parameter.endswith('_fname'):
                if value:
                    name = value.name
                    # path is created here to avoid creating empty folders
                    # for those modules without '_fname' parameters
                    pf = Path(data_dir, end_path)
                    pf.mkdir(exist_ok=True)
                    shutil.copy(value, Path(pf, name))
                    _p = Path(rel_data_dir, end_path, name)
                    new_mp[module][parameter] = _p

    return new_mp


def clean_rundir_according_to_restart(run_dir, restart_from=None):
    """
    Clean run directory according to restart parameter.

    Parameters
    ----------
    restart_from : None or int
        The module on which to restart the run. Discards all modules
        after this one (inclusive).
    """
    if restart_from is None:
        # prepares the run folders
        _p = Path(run_dir)
        if _p.exists() and len(list(_p.iterdir())) > 0:
            log.info(
                f"The {RUNDIR!r} {str(_p)!r} exists and is not empty. "
                "We can't work on it unless you provide the `--restart` "
                "option. If you want to start a run from scratch, "
                "indicate a new folder, or manually delete this one first, "
                "or use `--restart 0`."
                )
            sys.exit(get_goodbye_help())

    else:
        remove_folders_after_number(run_dir, restart_from)


def identify_modules(params):
    """Identify keys (headings) belogging to HADDOCK3 modules."""
    modules_keys = [
        k
        for k in params.keys()
        if get_module_name(k) in modules_category
        ]
    return modules_keys


def inject_in_modules(modules_params, key, value):
    """Inject a parameter in each module."""
    for params in modules_params.values():
        if key in params:
            raise ValueError(
                "key {key!r} already in {module!r} parameters. "
                "Can't inject."
                )
        params[key] = value


def validate_module_names_are_not_mispelled(params):
    """Validate headers are not mispelled."""
    module_names = sorted(modules_category.keys())
    for param_name, value in params.items():
        if isinstance(value, dict):
            if get_module_name(param_name) not in module_names:
                emsg = (
                    f"Module {param_name!r} is not a valid module name. "
                    f"Valid modules are: {', '.join(module_names)}."
                    )
                raise ValueError(emsg)


@with_config_error
def check_specific_validations(params):
    """Make specific validations."""
    # double check though this is confirmed already in the config reader
    v_rundir(params[RUNDIR])


# reading parameter blocks
def get_blocks(config):
    """
    Get parameter blocks.

    Block parameters follow the rule <preffix>_<something>_<digit>

    - part1_something1_1
    - part1_something2_1
    - part1_something_else_1
    - part1_something4_1

    When used to read the modules' default configuration we expect the
    <digit> to be only "_1". But having <digit> allows to identify
    blocks from the user configuration.

    We want to know:

    - part1
    - how many somethings exist (this defined the size of the block)
    - _1, defines the number of the block

    Returns
    -------
    dictionary
        In the form of:
        {("part1", "1"): {
            "something1", "something2", "something_else",
            "something4"}
    """
    splitted = (parameter.split("_") for parameter in config)
    parts = (
        _parts
        for _parts in splitted
        if len(_parts) > 2 and _parts[-1].isdigit()
        )

    blocks = {}
    for p in parts:
        new = blocks.setdefault((p[0], p[-1]), {})
        new.setdefault("counts", 0)
        new["counts"] += 1
        new.setdefault("mid", set())
        new["mid"].add("_".join(p[1:-1]))

    final_blocks = {k: v["mid"] for k, v in blocks.items() if v["counts"] > 1}

    return final_blocks


def read_blocks(eblocks, params):
    """
    Read the blocks from the user configuration file.

    Compare them with those defined in the default config.

    Parameters
    ----------
    eblocks : dict
        Parameter blocks identified in the `default.cfg` file using the
        :func:`get_blocks`.

    params : dict
        The user configuration dictionary.

    Returns
    -------
    ser
        A set of the new parameters that are allowed considering the
        groups found in the `default.cfg`.
    """
    pblocks = get_blocks(params)
    eblocks = {k: v for k, v in eblocks.items() if k[1] == "1"}

    enames = [_[0] for _ in eblocks]
    new = set()

    for block in pblocks:

        if block[0] not in enames:
            emsg = (
                f"The parameter block '{block[0]}_*_{block[1]}' "
                "is not a valid expandable parameter.")
            raise ConfigurationError(emsg)

        diff = pblocks[block].difference(eblocks[(block[0], "1")])
        if diff:
            emsg = (
                "These parameters do not belong to the block "
                f"'{block[0]}_*_{block[1]}': {', '.join(diff)}."
                )
            raise ConfigurationError(emsg)

        num_found = len(pblocks[block])  # len of the set of elements
        num_expected = len(eblocks[(block[0], "1")])

        if num_found < num_expected:
            emsg = (
                f"The parameter block '{block[0]}_*_{block[1]}' expects "
                f"{num_expected} parameters, but only {num_found} are present "
                "in the configuration file."
                )
            raise ConfigurationError(emsg)

        if num_found > num_expected:
            emsg = (
                f"The parameter block {block!r} expects {num_expected} "
                f"parameters, but {num_found} are present in the configuration "
                "file."
                )
            raise ConfigurationError(emsg)

        for p in params:
            try:
                pname, *_, pidx = p.split("_")
            except ValueError:
                continue
            if pname == block[0] and pidx == block[-1]:
                new.add(p)

    return new
