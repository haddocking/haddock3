"""Logic pertraining to preparing the run files and folders."""
import difflib
import importlib
import itertools as it
import os
import shutil
import string
import sys
import tarfile
from contextlib import contextmanager, suppress
from copy import copy
from functools import lru_cache, wraps
from pathlib import Path

from haddock import contact_us, haddock3_source_path, log
from haddock.core.defaults import RUNDIR, max_molecules_allowed
from haddock.core.exceptions import ConfigurationError, ModuleError
from haddock.gear.clean_steps import (
    UNPACK_FOLDERS,
    unpack_compressed_and_archived_files,
    update_unpacked_names,
    )
from haddock.gear.config import get_module_name
from haddock.gear.config import load as read_config
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
from haddock.gear.extend_run import (
    read_num_molecules_from_folder,
    renum_step_folders,
    )
from haddock.gear.greetings import get_goodbye_help
from haddock.gear.parameters import (
    config_mandatory_general_parameters,
    config_optional_general_parameters,
    config_optional_general_parameters_dict,
    )
from haddock.gear.preprocessing import process_pdbs, read_additional_residues
from haddock.gear.restart_run import remove_folders_after_number
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
    modules_names,
    non_mandatory_general_parameters_defaults,
    )
from haddock.modules.analysis import (
    confirm_resdic_chainid_length,
    modules_using_resdic,
    )


ALL_POSSIBLE_GENERAL_PARAMETERS = set.union(
    set(config_mandatory_general_parameters),
    set(non_mandatory_general_parameters_defaults),
    config_optional_general_parameters,
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
        extend_run=None,
        ):
    """
    Set up an HADDOCK3 run.

    This function sets up a HADDOCK3 considering the options `--restart`
    and `--extend-run`. The list of actions presented below does
    not necessary represents the exact order in which it happens.

    Always performed:

    #. read the user configuration file
    #. completes the user configuration file with the default values
       for the non-specified parameters
    #. validate the config file
       * confirm modules' names are correctly spelled
       * check if requested modules are installed
       * check additional validations
    #. validate modules' parameters
    #. copy input files to data/ directory
       * for ``--restart`` copies only after the restart number

    Performed when ``--restart``:

    #. remove folders after --restart number
    #. remove also folders from `data/` dir after the ``--restart`` num
    #. renumber step folders according to the number of modules

    Performed when ``--extend-run``:

    #. renumber step folders according to the number of modules

    Performed when start from scratch:

    #. check mandatory arguments are present in the config file
    #. check run-dir exists
    #. copy molecules to topology key (also in ``--restart``)
    #. populate topology parameters (also in ``--restart``)
    #. copy molecules to data dir

    Parameters
    ----------
    workflow_path : str or pathlib.Path
        The path to the configuration file.

    restart_from : int
        The step to restart the run from (inclusive).
        Defaults to None, which ignores this option.

    extend_run : str or Path
        The path created with `haddock3-copy` to start the run from.
        Defaults to None, which ignores this option.

    Returns
    -------
    tuple of two dicts
        A dictionary with the parameters for the haddock3 modules.
        A dictionary with the general run parameters.
    """
    # read the user config file from path
    params = read_config(workflow_path)

    # update default non-mandatory parameters with user params
    params = recursive_dict_update(
        config_optional_general_parameters_dict,
        params)

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

    # --extend-run configs do not define the run directory
    # in the config file. So we take it from the argument.
    if not_none(extend_run):
        with suppress(TypeError):
            extend_run = Path(extend_run)

        general_params[RUNDIR] = extend_run

    check_if_modules_are_installed(modules_params)
    check_specific_validations(general_params)

    # define starting conditions
    # consider a deeper refactor if additional conditions are implemented
    # @joaomcteixeira, 09 May 2022
    from_scratch = restart_from is None and extend_run is None
    scratch_rest0 = from_scratch or restart_from == 0
    restarting_from = not_none(restart_from)
    starting_from_copy = not_none(extend_run)

    if from_scratch:
        check_run_dir_exists(general_params[RUNDIR])

    if scratch_rest0:
        check_mandatory_argments_are_present(general_params)

    if restarting_from:
        remove_folders_after_number(general_params[RUNDIR], restart_from)
        _data_dir = Path(general_params[RUNDIR], "data")
        remove_folders_after_number(_data_dir, restart_from)

    if restarting_from or starting_from_copy:
        # get run files in folder
        step_folders = get_module_steps_folders(general_params[RUNDIR])

        log.info(
            'Uncompressing previous output files for folders: '
            f'{", ".join(step_folders)}'
            )
        # unpack the possible compressed and archived files
        _step_folders = (Path(general_params[RUNDIR], p) for p in step_folders)
        unpack_compressed_and_archived_files(
            _step_folders,
            general_params["ncores"],
            dec_all=True,
            )

    if starting_from_copy:
        num_steps = len(step_folders)
        _num_modules = len(modules_params)
        # has to consider the folders already present, plus the new folders
        # in the configuration file
        zero_fill.set_zerofill_number(num_steps + _num_modules)

        max_mols = read_num_molecules_from_folder(extend_run)

    else:
        copy_molecules_to_topology(
            general_params['molecules'],
            modules_params['topoaa.1'],
            )

        if len(modules_params["topoaa.1"]["molecules"]) > max_molecules_allowed:
            raise ConfigurationError("Too many molecules defined, max is {max_molecules_allowed}.")  # noqa: E501

        zero_fill.read(modules_params)

        populate_topology_molecule_params(modules_params["topoaa.1"])
        populate_mol_parameters(modules_params)

        max_mols = len(modules_params["topoaa.1"]["molecules"])

    if not from_scratch:
        _prev, _new = renum_step_folders(general_params[RUNDIR])
        renum_step_folders(Path(general_params[RUNDIR], "data"))
        if UNPACK_FOLDERS:  # only if there was any folder unpacked
            update_unpacked_names(_prev, _new, UNPACK_FOLDERS)
        update_step_contents_to_step_names(
            _prev,
            _new,
            general_params[RUNDIR],
            )

    validate_modules_params(modules_params, max_mols)

    # create datadir
    data_dir = create_data_dir(general_params[RUNDIR])

    if scratch_rest0:
        copy_molecules_to_data_dir(
            data_dir,
            modules_params["topoaa.1"],
            preprocess=general_params["preprocess"],
            )

    if starting_from_copy:
        copy_input_files_to_data_dir(data_dir, modules_params, start=num_steps)

    elif restarting_from:
        # copies only the input molecules needed
        _keys = list(modules_params.keys())
        _partial_params = {k: modules_params[k] for k in _keys[restart_from:]}
        copy_input_files_to_data_dir(
            data_dir,
            _partial_params,
            start=restart_from,
            )

    else:
        # copies everything
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
        module_name = get_module_name(module_name)
        defaults = _read_defaults(module_name)
        if not defaults:
            continue

        if module_name in modules_using_resdic:
            confirm_resdic_chainid_length(args)

        expandable_params = get_expandable_parameters(
            args,
            defaults,
            module_name,
            max_mols,
            )

        all_parameters = set.union(
            set(extract_keys_recursive(defaults)),
            set(non_mandatory_general_parameters_defaults.keys()),
            expandable_params,
            )

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
def copy_molecules_to_topology(molecules, topoaa_params):
    """Copy molecules to mandatory topology module."""
    topoaa_params['molecules'] = list(map(Path, transform_to_list(molecules)))


def copy_molecules_to_data_dir(data_dir, topoaa_params, preprocess=True):
    """
    Copy molecules to data directory and to topoaa parameters.

    Parameters
    ----------
    data_dir : Path
        The data/ directory inside the run directory. Must contain
        reference to the run directory.

    topoaa_params : dict
        A dictionary containing the topoaa parameters.

    preprocess : bool
        Whether to preprocess input molecules. Defaults to ``True``.
        See :py:mod:`haddock.gear.preprocessing`.
    """
    topoaa_dir = zero_fill.fill('topoaa', 0)

    # define paths
    data_topoaa_dir = Path(data_dir, topoaa_dir)
    data_topoaa_dir.mkdir(parents=True, exist_ok=True)
    rel_data_topoaa_dir = Path(data_dir.name, topoaa_dir)
    original_mol_dir = Path(data_dir, "original_molecules")

    new_molecules = []
    for molecule in copy(topoaa_params['molecules']):
        check_if_path_exists(molecule)

        mol_name = Path(molecule).name

        if preprocess:  # preprocess PDB files

            top_fname = topoaa_params.get("ligand_top_fname", False)
            new_residues = \
                read_additional_residues(top_fname) if top_fname else None

            new_pdbs = \
                process_pdbs(molecule, user_supported_residues=new_residues)

            # copy the original molecule
            original_mol_dir.mkdir(parents=True, exist_ok=True)
            original_mol = Path(original_mol_dir, mol_name)
            shutil.copy(molecule, original_mol)

            # write the new processed molecule
            new_pdb = os.linesep.join(new_pdbs[0])
            Path(data_topoaa_dir, mol_name).write_text(new_pdb)

        else:
            shutil.copy(molecule, Path(data_topoaa_dir, mol_name))

        new_molecules.append(Path(rel_data_topoaa_dir, mol_name))

    topoaa_params['molecules'] = copy(new_molecules)


def copy_input_files_to_data_dir(data_dir, modules_params, start=0):
    """
    Copy input files to data directory.

    Parameters
    ----------
    data_dir : Path
        The data/ directory inside the run directory. Must contain
        reference to the run directory.

    modules_params : dict
        A dictionary with the parameters of the modules. The paths to
        data in the dictionary are updated to the new paths copied
        to the data/ folder.

    start : int, default to 0
        The starting number of the step folders prefix.
    """
    rel_data_dir = data_dir.name
    for i, (module, params) in enumerate(modules_params.items(), start=start):
        end_path = Path(zero_fill.fill(get_module_name(module), i))
        for parameter, value in params.items():
            if parameter.endswith('_fname'):
                if value:
                    name = Path(value).name
                    # path is created here to avoid creating empty folders
                    # for those modules without '_fname' parameters
                    pf = Path(data_dir, end_path)
                    pf.mkdir(exist_ok=True)
                    check_if_path_exists(value)
                    target_path = Path(pf, name)
                    shutil.copy(value, target_path)
                    _p = Path(rel_data_dir, end_path, name)
                    modules_params[module][parameter] = _p
                    # account for input .tgz files
                    if name.endswith("tgz"):
                        log.info(f"Uncompressing tar {value}")
                        with tarfile.open(target_path) as fin:
                            fin.extractall(pf)


def check_run_dir_exists(run_dir):
    """Check whether the run directory exists."""
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


def validate_module_names_are_not_misspelled(params):
    """
    Validate module names are not misspelled in step definitions.

    Parameters
    ----------
    params : dict
        The user configuration file.
    """
    params_to_check = [
        get_module_name(param)
        for param, value in params.items()
        if isinstance(value, dict)
        ]

    validate_parameters_are_not_misspelled(
        params_to_check,
        reference_parameters=modules_names,
        )

    return


def validate_parameters_are_not_misspelled(params, reference_parameters):
    """Validate general parameters are not misspelled."""
    for param_name in params:
        if param_name not in reference_parameters:
            matched = fuzzy_match([param_name], reference_parameters)
            emsg = (
                f"Parameter {param_name!r} is not a valid general parameter,"
                f" did you mean {matched[0][1]!r}?"
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
    if get_module_name(module_name) == "topoaa":
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
        type_3 = type_simplest_ep[get_module_name(module_name)]
        allowed_params.update(read_simplest_expandable(type_3, user_config))

    _ = read_mol_parameters(user_config, type_4, max_mols=max_mols)
    allowed_params.update(_)

    return allowed_params


def populate_topology_molecule_params(topoaa):
    """Populate topoaa `molX` subdictionaries."""
    topoaa_dft = _read_defaults("topoaa.1")

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
        if not (mol in topoaa and "prot_segid" in topoaa[mol]):
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
    num_mols = range(1, len(modules_params["topoaa.1"]["molecules"]) + 1)
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


def update_step_contents_to_step_names(prev_names, new_names, folder):
    """
    Update step folder names in files after the `--restart` option.

    Runs over the folders defined in `new_names`.

    Parameters
    ----------
    prev_names : list
        List of step names to find in file contents.

    new_names : list
        List of new step names to replace `prev_names`. Both lists need
        to be synchronized. That is, the first index of `prev_names` should
        correspond to the old names of `new_names`.

    folder : str or Path
        Folder where the step folders are. Usually run directory or
        data directory.

    Returns
    -------
    None
        Save files in place.
    """
    for new_step in new_names:
        new_step_p = Path(folder, new_step)
        for file_ in new_step_p.iterdir():

            # goes recursive into the next folder
            if file_.is_dir():
                update_step_names_in_subfolders(file_, prev_names, new_names)

            else:
                update_step_names_in_file(file_, prev_names, new_names)


def update_step_names_in_subfolders(folder, prev_names, new_names):
    """
    Update step names in subfolders.

    Some modules may generate subfolders. This function update
    its files accordingly to the `--restart` feature.
    """
    for file_ in folder.iterdir():
        if file_.is_dir():
            update_step_names_in_subfolders(file_, prev_names, new_names)
        else:
            update_step_names_in_file(file_, prev_names, new_names)
    return


def update_step_names_in_file(file_, prev_names, new_names):
    """Update step names in file following the `--restart` option."""
    text = file_.read_text()
    for s1, s2 in zip(prev_names, new_names):
        text = text.replace(s1, s2)
    file_.write_text(text)
    return
