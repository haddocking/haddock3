"""Logic pertraining to preparing the run files and folders."""
import difflib
import importlib
import itertools as it
import os
import shutil
import string
import sys
from contextlib import contextmanager, suppress
from functools import lru_cache, wraps
from pathlib import Path

from haddock import contact_us, haddock3_source_path, log
from haddock.core.defaults import RUNDIR, max_molecules_allowed
from haddock.core.exceptions import ConfigurationError, ModuleError
from haddock.gear.config_reader import get_module_name, read_config
from haddock.gear.expandable_parameters import (
    get_mol_parameters,
    get_multiple_index_groups,
    get_single_index_groups,
    is_mol_parameter,
    read_mol_parameters,
    read_multiple_idx_groups_user_config,
    read_simplest_expandable,
    read_single_idx_groups_user_config,
    remove_trail_idx,
    type_simplest_ep,
    )
from haddock.gear.greetings import get_goodbye_help
from haddock.gear.parameters import config_mandatory_general_parameters
from haddock.gear.restart_run import remove_folders_after_number
from haddock.gear.restart_from_copy import read_num_molecules_from_folder
from haddock.gear.validations import v_rundir
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.gear.zerofill import zero_fill
from haddock.libs.libfunc import not_none
from haddock.libs.libutil import (
    extract_keys_recursive,
    recursive_dict_update,
    remove_dict_keys,
    transform_to_list,
    )
from haddock.modules import (
    get_module_steps_folders,
    modules_category,
    non_mandatory_general_parameters_defaults,
    )
from haddock.modules.analysis import (
    confirm_resdic_chainid_length,
    modules_using_resdic,
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


@lru_cache
def _read_defaults(module_name):
    """Read the defaults.yaml given a module name."""
    module_name_ = get_module_name(module_name)
    pdef = Path(
        haddock3_source_path,
        'modules',
        modules_category[module_name_],
        module_name_,
        'defaults.yaml',
        ).resolve()

    return read_from_yaml_config(pdef)


def setup_run(
        workflow_path,
        restart_from=None,
        restart_from_copy=None,
        ):
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

    with suppress(TypeError):
        restart_from_copy = Path(restart_from_copy)

    if not_none(restart_from_copy):
        params[RUNDIR] = restart_from_copy

    if restart_from_copy is None:
        check_mandatory_argments_are_present(params)

    validate_module_names_are_not_mispelled(params)
    check_specific_validations(params)

    # update default non-mandatory parameters with user params
    params = recursive_dict_update(
        non_mandatory_general_parameters_defaults,
        params)

    if restart_from_copy is None:
        clean_rundir_according_to_restart(params[RUNDIR], restart_from)

    # copy molecules parameter to topology module
    if restart_from_copy is None:
        copy_molecules_to_topology(params)
        if len(params["topoaa"]["molecules"]) > max_molecules_allowed:
            raise ConfigurationError("Too many molecules defined, max is {max_molecules_allowed}.")  # noqa: E501

    # separate general from modules parameters
    _modules_keys = identify_modules(params)
    general_params = remove_dict_keys(params, _modules_keys)
    modules_params = remove_dict_keys(params, list(general_params.keys()))

    if not_none(restart_from_copy):
        num_steps = len(get_module_steps_folders(restart_from_copy))
        _num_modules = len(modules_params)
        # has to consider the folders already present, plus the new folders
        # in the configuration file
        zero_fill.set_zerofill_number(num_steps + _num_modules)
    else:
        zero_fill.read(modules_params)

    # populate topology molecules
    if restart_from_copy is None:
        populate_topology_molecule_params(modules_params["topoaa"])
        populate_mol_parameters(modules_params)

    # validations
    if restart_from_copy is None:
        max_mols = len(modules_params["topoaa"]["molecules"])
    else:
        max_mols = read_num_molecules_from_folder(restart_from_copy)

    validate_modules_params(modules_params, max_mols)
    check_if_modules_are_installed(modules_params)

    # create datadir
    data_dir = create_data_dir(general_params[RUNDIR])

    if restart_from_copy is None:
        copy_molecules_to_data_dir(data_dir, modules_params["topoaa"])

    if not_none(restart_from_copy):
        copy_input_files_to_data_dir(data_dir, modules_params, start=num_steps)
    else:
        copy_input_files_to_data_dir(data_dir, modules_params)

    # return the modules' parameters and general parameters separately
    return modules_params, general_params


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
def validate_modules_params(modules_params, max_mols):
    """
    Validate individual parameters for each module.

    Raises
    ------
    ConfigError
        If there is any parameter given by the user that is not defined
        in the defaults.cfg of the module.
    """
    for module_name, args in modules_params.items():
        defaults = _read_defaults(module_name)
        if not defaults:
            return

        if module_name in modules_using_resdic:
            confirm_resdic_chainid_length(args)

        expandable_params = get_expandable_parameters(
            args,
            defaults,
            module_name,
            max_mols,
            )

        all_parameters = \
            set.union(set(extract_keys_recursive(defaults)),
                      set(config_mandatory_general_parameters),
                      set(non_mandatory_general_parameters_defaults.keys()),
                      expandable_params)

        diff = set(extract_keys_recursive(args)) - all_parameters

        if diff:
            matched = fuzzy_match(diff, all_parameters)

            def pretty_print(match):
                return f" * \'{match[0]}\' did you mean \'{match[1]}\'?"

            _msg = (
                'The following parameters do not match any expected '
                f'parameters for module {module_name!r}: {os.linesep}'
                f'{os.linesep.join(map(pretty_print, matched))}.'
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


def copy_molecules_to_data_dir(data_dir, topoaa_params):
    """Copy files to data directory."""
    # this line must be synchronized with create_data_dir()
    rel_data_dir = data_dir.name

    topoaa_dir = zero_fill.fill('topoaa', 0)
    for i, molecule in enumerate(topoaa_params['molecules']):
        end_path = Path(data_dir, topoaa_dir)
        end_path.mkdir(parents=True, exist_ok=True)
        name = Path(molecule).name
        check_if_path_exists(molecule)
        shutil.copy(molecule, Path(end_path, name))
        topoaa_params['molecules'][i] = Path(rel_data_dir, topoaa_dir, name)


def copy_input_files_to_data_dir(data_dir, modules_params, start=0):
    # topology always starts with 0
    rel_data_dir = data_dir.name
    for i, (module, params) in enumerate(modules_params.items(), start=start):
        end_path = Path(zero_fill.fill(get_module_name(module), i))
        for parameter, value in params.items():
            if parameter.endswith('_fname'):
                if value:
                    name = value.name
                    # path is created here to avoid creating empty folders
                    # for those modules without '_fname' parameters
                    pf = Path(data_dir, end_path)
                    pf.mkdir(exist_ok=True)
                    check_if_path_exists(value)
                    shutil.copy(value, Path(pf, name))
                    _p = Path(rel_data_dir, end_path, name)
                    modules_params[module][parameter] = _p


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
    """Identify keys (headings) belonging to HADDOCK3 modules."""
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
    """Validate headers are not misspelled."""
    module_names = sorted(modules_category.keys())
    for param_name, value in params.items():
        if isinstance(value, dict):
            module_name = get_module_name(param_name)
            if module_name not in module_names:
                matched = fuzzy_match([module_name], module_names)
                emsg = (
                    f"Module {param_name!r} is not a valid module name,"
                    f" did you mean {matched[0][1]}?. "
                    f"Valid modules are: {', '.join(module_names)}."
                    )
                raise ValueError(emsg)


@with_config_error
def check_specific_validations(params):
    """Make specific validations."""
    # double check though this is confirmed already in the config reader
    v_rundir(params[RUNDIR])


def get_expandable_parameters(user_config, defaults, module_name, max_mols):
    """
    Get configuration expandable blocks.

    Parameters
    ----------
    user_config : dict
        The user configuration file for a module.

    defaults : dict
        The default configuration file defined for the module.

    module_name : str
        The name the module being processed.

    max_mols : int
        The max number of molecules allowed.
    """
    # the topoaa module is an exception because it has subdictionaries
    # for the `mol` parameter. Instead of defining a general recursive
    # function, I decided to add a simple if/else exception.
    # no other module should have subdictionaries has parameters
    if module_name == "topoaa":
        ap = set()  # allowed_parameters
        ap.update(_get_expandable(user_config, defaults, module_name, max_mols))
        for i in range(1, max_mols + 1):
            key = f"mol{i}"
            with suppress(KeyError):
                ap.update(
                    _get_expandable(
                        user_config[key],
                        defaults["mol1"],
                        module_name,
                        max_mols,
                        )
                    )

        return ap

    elif module_name in modules_using_resdic:
        ep = _get_expandable(user_config, defaults, module_name, max_mols)
        for _param in user_config.keys():
            if _param.startswith("resdic_"):
                ep.add(_param)
        return ep

    else:
        return _get_expandable(user_config, defaults, module_name, max_mols)


# reading parameter blocks
def _get_expandable(user_config, defaults, module_name, max_mols):
    type_1 = get_single_index_groups(defaults)
    type_2 = get_multiple_index_groups(defaults)
    type_4 = get_mol_parameters(defaults)

    allowed_params = set()
    allowed_params.update(read_single_idx_groups_user_config(user_config, type_1))  # noqa: E501
    allowed_params.update(read_multiple_idx_groups_user_config(user_config, type_2))  # noqa: E501

    with suppress(KeyError):
        type_3 = type_simplest_ep[module_name]
        allowed_params.update(read_simplest_expandable(type_3, user_config))

    _ = read_mol_parameters(user_config, type_4, max_mols=max_mols)
    allowed_params.update(_)

    return allowed_params


def populate_topology_molecule_params(topoaa):
    """Populate topoaa `molX` subdictionaries."""
    topoaa_dft = _read_defaults("topoaa")

    # list of possible prot_segids
    uppers = list(string.ascii_uppercase)[::-1]

    # removes from the list those prot_segids that are already defined
    for param in topoaa:
        if param.startswith("mol") and param[3:].isdigit():
            with suppress(KeyError):
                uppers.remove(topoaa[param]["prot_segid"])

    # populates the prot_segids just for those that were not defined
    # in the user configuration file. Other parameters are populated as
    # well. `prot_segid` is the only one differing per molecule.
    for i in range(1, len(topoaa["molecules"]) + 1):
        mol = f"mol{i}"
        if not(mol in topoaa and "prot_segid" in topoaa[mol]):
            topoaa_dft["mol1"]["prot_segid"] = uppers.pop()

        topoaa[mol] = recursive_dict_update(
            topoaa_dft["mol1"],
            topoaa[mol] if mol in topoaa else {},
            )
    return


def populate_mol_parameters(modules_params):
    """
    Populate modules subdictionaries with the needed molecule `mol_` parameters.

    The `mol_` prefixed parameters is a subclass of the expandable parameters.

    See `gear.expandable_parameters`.

    Modules require these parameters to be repeated for the number of input
    molecules.

    This function adds `mol_` parameters to the user input parameters,
    one per each `molecule`.

    Parameters
    ----------
    modules_params : dict
        A dictionary containing only modules' keys:subdictionaries
        parameters. That is, without the general parameters.

    Returns
    -------
    None
        Alter the dictionary in place.
    """
    # the starting number of the `mol_` parameters is 1 by CNS definition.
    num_mols = range(1, len(modules_params["topoaa"]["molecules"]) + 1)
    for module_name, _ in modules_params.items():

        # read the modules default parameters
        defaults = _read_defaults(module_name)

        # if there are no `mol_` parameters in the modules default values,
        # the `mol_params` generator will be empty and the for-loop below
        # won't run.
        mol_params = (p for p in list(defaults.keys()) if is_mol_parameter(p))

        for param, i in it.product(mol_params, num_mols):
            param_name = remove_trail_idx(param)

            # the `setdefault` grants that the value is only added if
            # the parameter is not present.
            modules_params[module_name].setdefault(
                f"{param_name}_{i}",
                defaults[param],
                )
    return


def check_if_path_exists(path):
    """
    Check if a path exists and raises an error if it does not exist.

    For example given this path "../config/project_01/file.txt" it would find
    the following path "../config/project-01".

    Parameters
    ----------
    path : AnyStr | PathLike
        The path to check.

    Returns
    -------
    None
        If the path does exist.

    Raises
    ------
    ValueError
        If the path does not exist.
    """
    path = os.path.normpath(path)
    if os.path.exists(path):
        return None

    reconstituted_path = "./"
    error = ("", "", "")
    elements = Path(path).parts
    if elements[0] == ".":
        elements = elements[1:]
    for part in elements:
        next_folder = Path(reconstituted_path, part)
        if not next_folder.exists():
            error = (reconstituted_path, fuzzy_match([part],
                     os.listdir(reconstituted_path))[0][1], part)
            break
        reconstituted_path = next_folder

    msg = (f"The following file could not be found: \'{path}\'. "
           f"In the folder \'{error[0]}\' the following \'{error[1]}\' "
           f"is the closest match to the supplied \'{error[2]}\', did "
           "you mean to open this?")
    raise ValueError(msg)


def fuzzy_match(user_input, possibilities):
    """
    Find the closest possibility to the user supplied input.

    Parameters
    ----------
    user_input : list(string)
        List of strings with the faulty input given by the user.
    possibilities : list(string)
        List of strings with all possible options that would be
        valid in this context.

    Returns
    -------
    list(string, string)
        The closest string from the possibilities to each string of the
        `user_input`. With as first element of the tuple the user_input
        string, and as second element the matched possibility.
    """
    results = list()

    for user_word in transform_to_list(user_input):
        best = (-1, "")
        for possibility in possibilities:
            distance = difflib.SequenceMatcher(a=user_word, b=possibility).ratio()  # noqa: E501
            if distance > best[0]:
                best = (distance, possibility)
        results += [(user_word, best[1])]

    return results
