"""Logic pertraining to preparing the run files and folders."""
import difflib
import importlib
import itertools as it
import json
import os
import re
import shutil
import sys
import tarfile
from contextlib import contextmanager, suppress
from copy import copy, deepcopy
from functools import lru_cache, wraps
from pathlib import Path, PosixPath

from haddock import EmptyPath, contact_us, haddock3_source_path, log
from haddock.core.defaults import RUNDIR, max_molecules_allowed, DATA_DIRNAME
from haddock.core.exceptions import (
    ConfigurationError,
    ModuleError,
    DependencyError,
    )
from haddock.core.typing import (
    Any,
    Callable,
    FilePath,
    Generator,
    Iterable,
    Optional,
    ParamDict,
    ParamMap,
    Union,
)
from haddock.gear.clean_steps import (
    UNPACK_FOLDERS,
    unpack_compressed_and_archived_files,
    update_unpacked_names,
)
from haddock.gear.config import get_module_name
from haddock.gear.config import load as read_config
from haddock.gear.config import save as save_config
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
from haddock.libs.libio import make_writeable_recursive
from haddock.libs.libutil import (
    extract_keys_recursive,
    recursive_convert_paths_to_strings,
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

# Dict mapping string types (in default.yaml) into python3 objects types
TYPES_MAPPER = {
    "boolean": (bool,),
    "bool": (bool,),
    "integer": (
        int,
        float,
    ),
    "int": (
        int,
        float,
    ),
    "long": (
        float,
        int,
    ),
    "float": (
        float,
        int,
    ),
    "double": (float,),
    "string": (str,),
    "str": (str,),
    "list": (list,),
    "array": (list,),
    "path": (
        str,
        Path,
        PosixPath,
        EmptyPath,
    ),
    "file": (
        str,
        Path,
        PosixPath,
        EmptyPath,
    ),
}


@contextmanager
def config_key_error() -> Generator[None, None, None]:
    """Raise ConfigurationError on KeyError."""
    try:
        yield
    except KeyError as err:
        msg = f"Expected {err.args[0]!r} parameter in configuration file."
        raise ConfigurationError(msg) from err


def with_config_error(func: Callable[..., Any]) -> Callable[..., Any]:
    """Add config error context."""

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:
        with config_key_error():
            return func(*args, **kwargs)

    return wrapper


@lru_cache
def _read_defaults(module_name, default_only=True):
    """Read the defaults.yaml given a module name.

    Parameters
    ----------
    module_name : str
        Name of the HADDOCK3 module
    default_only : bool
        If True, only the default value of each parameters found in the
        default.yaml will be return, else all information are retrieved.

    Return
    ------
    mod_default_config : dict
        A dict holding default value for each parameters of the module if
        default_only == True,
        Else a dict of dict, contraining all information present in the
        default.yaml file of the module.
    """
    module_name_ = get_module_name(module_name)
    pdef = Path(
        haddock3_source_path,
        "modules",
        modules_category[module_name_],
        module_name_,
        "defaults.yaml",
    ).resolve()
    mod_default_config = read_from_yaml_config(pdef, default_only=default_only)
    return mod_default_config


def setup_run(
    workflow_path: FilePath,
    restart_from: Optional[int] = None,
    extend_run: Optional[FilePath] = None,
) -> tuple[ParamDict, ParamDict]:
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
    config_files = read_config(workflow_path)

    # update default non-mandatory parameters with user params
    params = recursive_dict_update(
        config_optional_general_parameters_dict, config_files["final_cfg"]
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
        check_CNS_usage(modules_params)

    if scratch_rest0:
        check_mandatory_argments_are_present(general_params)

    if restarting_from:
        remove_folders_after_number(general_params[RUNDIR], restart_from)
        _data_dir = Path(general_params[RUNDIR], DATA_DIRNAME)
        remove_folders_after_number(_data_dir, restart_from)

    if restarting_from or starting_from_copy:
        # get run files in folder
        step_folders = get_module_steps_folders(general_params[RUNDIR])

        log.info(
            "Uncompressing previous output files for folders: "
            f'{", ".join(step_folders)}'
        )
        # unpack the possible compressed and archived files
        _step_folders = (Path(general_params[RUNDIR], p) for p in step_folders)
        unpack_compressed_and_archived_files(
            _step_folders,
            general_params["ncores"],
            dec_all=True,
        )

    first_module_id = list(modules_params.keys())[0]
    # Here we check if topoaa is the first module in the workflow.
    # If it is, we gather the parameters of topoaa in the topology_params,
    # as equired by the function populate_topology_molecule_params(),
    # to map the molX parameters info to input molecules.
    # Without this, we loose the information.
    if first_module_id == "topoaa.1":
        topology_params = modules_params[first_module_id]
    else:
        # If not, just fake an empty set of parameters
        topology_params = {}

    if starting_from_copy:
        num_steps = len(step_folders)
        _num_modules = len(modules_params)
        # has to consider the folders already present, plus the new folders
        # in the configuration file
        zero_fill.set_zerofill_number(num_steps + _num_modules)

        max_mols = read_num_molecules_from_folder(extend_run)

    else:
        copy_molecules_to_topology(
            general_params["molecules"],
            topology_params,
        )

        max_mols = len(topology_params["molecules"])
        if max_mols > max_molecules_allowed:
            raise ConfigurationError(
                f"Too many molecules defined, max is {max_molecules_allowed}."
                )

        zero_fill.read(modules_params)

        populate_topology_molecule_params(topology_params)
        populate_mol_parameters(modules_params, topology_params)

    if not from_scratch:
        _prev, _new = renum_step_folders(general_params[RUNDIR])
        renum_step_folders(Path(general_params[RUNDIR], DATA_DIRNAME))
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

    # Add workflow configuration file in data directory
    enhanced_haddock_params = deepcopy(general_params)
    enhanced_haddock_params.update(modules_params)
    config_files["enhanced_haddock_params"] = enhanced_haddock_params
    config_saves = save_configuration_files(config_files, data_dir)  # noqa : F841

    if scratch_rest0:
        copy_molecules_to_data_dir(
            data_dir,
            topology_params,
            first_module_id,
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

    # grant write permissions to data/ dir
    make_writeable_recursive(data_dir)

    # return the modules' parameters and general parameters separately
    return modules_params, general_params


def save_configuration_files(configs: dict, datadir: Union[str, Path]) -> dict:
    """Write a copy of configuration files (GitHub issue #578).

    Parameters
    ----------
    configs : dict
        Dictionnary holding the various configuration files
        ['raw_input, 'cleaned_input', 'loaded_cleaned_input',
        'final_cfg', 'enhanced_haddock_params']
    datadir : str or :py:class:`libpath.Path`
        Directory where to write the configuration.

    Return
    ------
    added_files : dict
        Dictionary of paths leading to saved configuration files.
    """
    # Create directory
    confpaths = Path(datadir, "configurations/")
    confpaths.mkdir(parents=True, exist_ok=True)
    # Initiate files data
    infofile = {
        "raw_input": (
            "An untouched copy of the raw input file, "
            "as provided by the user."
            ),
        "cleaned_input": (
            "Pre-parsed input file where (eventually) "
            "some indexing and modifications were "
            "applied to ensure further processing."
        ),
        "enhanced_haddock_params": (
            "Final input file with detailed default parameters."
            ),
        }
    added_files = {}
    # Set list of configurations that wish to be saved
    list_save_conf = [
        "raw_input",
        "cleaned_input",
        "enhanced_haddock_params",
    ]

    # Loop over configuration files
    for confname in list_save_conf:
        try:
            confdt = configs[confname]
        except KeyError:
            continue

        toml_fpath = Path(confpaths, f"{confname}.toml")
        if isinstance(confdt, dict):
            # Save toml version
            save_config(confdt, toml_fpath)
            # Save json version
            json_fpath = Path(confpaths, f"{confname}.json")
            jsonconf = recursive_convert_paths_to_strings(confdt)
            with open(json_fpath, "w", encoding="utf-8") as f:
                json.dump(jsonconf, f, indent=4)
        else:
            with open(toml_fpath, "w") as f:
                f.write(confdt)

        added_files[confname] = toml_fpath

    # Add README to help user
    readmepath = Path(confpaths, "README.txt")
    with open(readmepath, "w") as f:
        f.write(
            f"{'#'*80}\n# Information about configuration "
            f"files present in the same directory #\n{'#'*80}\n"
        )
        for confname in list_save_conf:
            if confname not in added_files.keys():
                continue
            f.write(f'"{added_files[confname]}": {infofile[confname]}\n')
    added_files["readme"] = readmepath

    # Return mapper to added files
    return added_files


def validate_params(params):
    """
    Validate the parameter file.

    #1 : checks for mandatory parameters
    #2 : checks for correct modules
    """
    check_mandatory_argments_are_present(params)
    validate_modules_names(params)


def check_mandatory_argments_are_present(params: Iterable[str]) -> None:
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
def validate_modules_names(params: Iterable[str]) -> None:
    """Validate all modules names are spelled correctly."""
    keys = (
        set(params)
        - set(config_mandatory_general_parameters)
        - set(non_mandatory_general_parameters_defaults)
    )

    for module in keys:
        if get_module_name(module) not in modules_category.keys():
            _msg = (
                f"Module {module} not found in HADDOCK3 library. "
                "Please refer to the list of available modules at: "
                "DOCUMENTATION-LINK"
            )
            raise ConfigurationError(_msg)


@with_config_error
def validate_modules_params(modules_params: ParamMap, max_mols: int) -> None:
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

            def pretty_print(match: tuple[str, str]) -> str:
                return f" * '{match[0]}' did you mean '{match[1]}'?"

            _msg = (
                "The following parameters do not match any expected "
                f"parameters for module {module_name!r}: {os.linesep}"
                f"{os.linesep.join(map(pretty_print, matched))}."
            )
            raise ConfigurationError(_msg)

        # Now that existence of the parameter was checked,
        # it is time to validate the type and value taken by this param
        validate_module_params_values(module_name, args)

        # Update parameters with defaults ones
        _missing_defaults = {
            param: default
            for param, default in defaults.items()
            if param not in args.keys()
        }
        args.update(_missing_defaults)


def validate_module_params_values(module_name: str, args: dict) -> None:
    """Validate individual parameters for a module.

    Parameters
    ----------
    module_name :
        Name of the module to be analyzed
    args : dict
        Dictionnary of key/value present in user config file for a module

    Raises
    ------
    ConfigError
        If there is any parameter given by the user that is not following the
        types or ranges/choices allowed in the defaults.cfg of the module.
    """
    # Load all parameters information from 'defaults.cfg'
    default_conf_params = _read_defaults(module_name, default_only=False)
    # Loop over user queried parameters keys/values
    for key, val in args.items():
        validate_value(default_conf_params, key, val)


def validate_value(default_yaml: dict, key: str, value: Any) -> None:
    """Validate queried value for a specific parameter of a module.

    Parameters
    ----------
    default_yaml : dict
        Dictionnary of key/value present in user config file for a module
    key : str
        Key to be analyzed
    value : [bool, int, float, str, list]
        The provided value for a parameter

    Raises
    ------
    ConfigError
        If there is any parameter given by the user that is not following the
        types or ranges/choices allowed in the defaults.cfg of the module.
    """
    # Special case for molecules...
    if key not in default_yaml.keys():
        return
    if "group" in default_yaml[key].keys():
        if default_yaml[key]["group"] == "molecules":
            return

    # Series of checks
    if _msg := validate_param_type(default_yaml[key], value):
        raise ConfigurationError(f'Config error for parameter "{key}": {_msg}')
    if _msg := validate_param_range(default_yaml[key], value):
        raise ConfigurationError(f'Config error for parameter "{key}": {_msg}')
    # FIXME: add more checks here ?


def validate_param_type(param: dict, val: Any) -> Optional[str]:
    """Check if provided parameter type is similar to defined defaults ones.

    Parameters
    ----------
    param : dict
        Dictionnary of key/value present in user config file for a module
    val : [bool, int, float, str, list]
        The provided value for a parameter

    Return
    ------
    May return a string explaining the issue with the provided value type
    """
    if "type" not in param.keys():
        return
    # Load type from string to python class
    try:
        allowed_types = TYPES_MAPPER[param["type"]]
    except KeyError as e:
        return f'Unrecognized default type "{e}"'
    # Check if both of the same type
    if (query_type := type(val)) not in allowed_types:
        return (
            f'Wrong provided type "{query_type}" with value "{val}". '
            f'It should be of type "{param["type"]}" '
            f'(e.g: {param["default"]})'
        )


def validate_param_range(param: dict, val: Any) -> Optional[str]:
    """Check if provided value is in range/choices defined for this parameter.

    Parameters
    ----------
    param : dict
        Dictionnary of key/value present in user config file for a module
    val : [bool, int, float, str, list]
        The provided value for a parameter

    Return
    ------
    May return a string explaining the issue with the provided value range
    """
    # Case for choices
    if "choices" in param.keys():
        if val not in (choices := param["choices"]):
            return f'Value "{val}" is not among the accepted choices: {choices}'  # noqa : E501
    # Case for ranges
    elif "min" in param.keys() and "max" in param.keys():
        if val < param["min"] or val > param["max"]:
            return (
                f'Value "{val}" is not in the allowed boundaries '
                f'ranging from {param["min"]} to {param["max"]}'
            )
    # Case for lists
    elif "minitems" in param.keys() and "maxitems" in param.keys():
        if (nbitems := len(val)) < param["minitems"] or nbitems > param["maxitems"]:
            _desc = "under" if nbitems < param["minitems"] else "exceeding"
            return (
                f'Number of items found in "{val}" is {_desc} the '
                "permitted limit. It should range from "
                f'{param["minitems"]} to {param["maxitems"]}'
            )
    # Case for strings
    elif "minchars" in param.keys() and "maxchars" in param.keys():
        if (nbchars := len(val)) < param["minchars"] or nbchars > param["maxchars"]:
            _desc = "under" if nbchars < param["minitems"] else "exceeding"
            return (
                f'Number of items found in "{val}" is {_desc} the '
                "permitted limit. It should range from "
                f'{param["minchars"]} to {param["maxchars"]}'
            )


def check_if_modules_are_installed(params: ParamMap) -> None:
    """Validate if third party-libraries are installed."""
    for module_name in params.keys():
        module_import_name = ".".join(
            [
                "haddock",
                "modules",
                modules_category[get_module_name(module_name)],
                get_module_name(module_name),
            ]
        )
        module_lib = importlib.import_module(module_import_name)
        try:
            module_lib.HaddockModule.confirm_installation()
        except Exception as err:
            _msg = (
                "A problem occurred while assessing if module "
                f"{module_name!r} is installed in your system. Have you "
                "installed the packages required to run this module? If "
                f"yes, write us at {contact_us!r} describing your system "
                "and the problems you are facing. If not, please install "
                "the required packages to use the module."
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
def create_data_dir(run_dir: FilePath) -> Path:
    """
    Create initial files for HADDOCK3 run.

    Returns
    -------
    pathlib.Path
        A path referring only to 'data'.
    """
    data_dir = Path(run_dir, DATA_DIRNAME)
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


@with_config_error
def copy_molecules_to_topology(
    molecules: Iterable[str], topoaa_params: ParamMap
) -> None:
    """Copy molecules to mandatory topology module."""
    topoaa_params["molecules"] = list(map(Path, transform_to_list(molecules)))


def copy_molecules_to_data_dir(
        data_dir: Path,
        topoaa_params: ParamMap,
        _first_module_name: str,
        preprocess: bool = True,
        ) -> None:
    """
    Copy molecules to data directory and to topoaa parameters.

    Parameters
    ----------
    data_dir : Path
        The data/ directory inside the run directory. Must contain
        reference to the run directory.

    topoaa_params : dict
        A dictionary containing the topoaa parameters.
    
    _first_module_name : str
        Name of the first module used in the workflow.

    preprocess : bool
        Whether to preprocess input molecules. Defaults to ``True``.
        See :py:mod:`haddock.gear.preprocessing`.
    """
    # Removes digit from module name
    # Build regex to capture '<name>.<digit>'
    name_digit_regex = re.compile(r"(\w+)(\.\d+)?")
    first_module_name: str = "input_molecules"
    if match := name_digit_regex.search(_first_module_name):
        first_module_name = match.group(1)
    topoaa_dir = zero_fill.fill(first_module_name, 0)

    # define paths
    data_topoaa_dir = Path(data_dir, topoaa_dir)
    data_topoaa_dir.mkdir(parents=True, exist_ok=True)
    rel_data_topoaa_dir = Path(data_dir.name, topoaa_dir)
    original_mol_dir = Path(data_dir, "original_molecules")

    # Init new molecule holder to be filled with relative paths
    new_molecules: list[Path] = []
    # Loop over input molecules
    for molecule in copy(topoaa_params["molecules"]):
        check_if_path_exists(molecule)
        mol_name = Path(molecule).name
        # preprocess PDB files
        if preprocess:
            # copy the un-processed molecule (for later checks)
            original_mol_dir.mkdir(parents=True, exist_ok=True)
            original_mol = Path(original_mol_dir, mol_name)
            shutil.copy(molecule, original_mol)
            # Gather potential user-provided topology file
            top_fname = topoaa_params.get("ligand_top_fname", False)
            new_residues = read_additional_residues(top_fname) if top_fname else None
            # Do the pre-processing of file
            new_pdbs = process_pdbs(molecule, user_supported_residues=new_residues)
            # write the new processed molecule
            new_pdb = os.linesep.join(new_pdbs[0])
            Path(data_topoaa_dir, mol_name).write_text(new_pdb)
        # Do not preprocess
        else:
            # Create a copy of input molecules into `data/0_firstmodule`
            shutil.copy(molecule, Path(data_topoaa_dir, mol_name))
        # Create relative path of the molecule
        data_dir_molecule_relpath = Path(rel_data_topoaa_dir, mol_name)
        new_molecules.append(data_dir_molecule_relpath)
    # Modify molecules parameters to point relative path of copied files
    topoaa_params["molecules"] = copy(new_molecules)


def copy_input_files_to_data_dir(
    data_dir: Path, modules_params: ParamMap, start: int = 0
) -> None:
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
            if parameter.endswith("_fname"):
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


def check_run_dir_exists(run_dir: FilePath) -> None:
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


def identify_modules(params: Iterable[str]) -> list[str]:
    """Identify keys (headings) belonging to HADDOCK3 modules."""
    modules_keys = [
        param_name for param_name in params
        if get_module_name(param_name) in modules_category
        ]
    return modules_keys


def inject_in_modules(modules_params: ParamMap, key: Any, value: Any) -> None:
    """Inject a parameter in each module."""
    for params in modules_params.values():
        if key in params:
            raise ValueError(
                "key {key!r} already in {module!r} parameters. " "Can't inject."
            )
        params[key] = value


def validate_module_names_are_not_misspelled(params: ParamMap) -> None:
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


def validate_parameters_are_not_misspelled(
    params: Iterable[str], reference_parameters: Iterable[str]
) -> None:
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
def check_specific_validations(params: ParamMap) -> None:
    """Make specific validations."""
    # double check though this is confirmed already in the config reader
    v_rundir(params[RUNDIR])


def get_expandable_parameters(
    user_config: ParamMap, defaults: ParamMap, module_name: str, max_mols: int
) -> set[str]:
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
    # no other module should have subdictionaries as parameters
    if get_module_name(module_name) == "topoaa":
        ap: set[str] = set()  # allowed_parameters
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
def _get_expandable(
    user_config: ParamMap, defaults: ParamMap, module_name: str, max_mols: int
) -> set[str]:
    type_1 = get_single_index_groups(defaults)
    type_2 = get_multiple_index_groups(defaults)
    type_4 = get_mol_parameters(defaults)

    allowed_params: set[str] = set()
    allowed_params.update(
        read_single_idx_groups_user_config(user_config, type_1)
    )  # noqa: E501
    allowed_params.update(
        read_multiple_idx_groups_user_config(user_config, type_2)
    )  # noqa: E501

    with suppress(KeyError):
        type_3 = type_simplest_ep[get_module_name(module_name)]
        allowed_params.update(read_simplest_expandable(type_3, user_config))

    _ = read_mol_parameters(user_config, type_4, max_mols=max_mols)
    allowed_params.update(_)

    return allowed_params


def populate_topology_molecule_params(topology_params: ParamMap) -> None:
    """Populate topoaa `molX` subdictionaries.

    Parameters
    ----------
    topology_params : ParamMap
        Dictionary of parameter with their values.
        Possibily parameters from topoaa module.
        If not, nothing will happen
    """
    topoaa_dft = _read_defaults("topoaa.1")

    for i in range(1, len(topology_params["molecules"]) + 1):
        mol = f"mol{i}"

        topology_params[mol] = recursive_dict_update(
            topoaa_dft["mol1"],
            topology_params[mol] if mol in topology_params else {},
        )
    return


def populate_mol_parameters(
        modules_params: ParamMap,
        topology_params: ParamMap,
        ) -> None:
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
    num_mols = range(1, len(topology_params["molecules"]) + 1)
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


def check_if_path_exists(path: FilePath) -> None:
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

    reconstituted_path: str = "./"
    error = ("", "", "")
    elements = Path(path).parts
    if elements[0] == ".":
        elements = elements[1:]
    for part in elements:
        next_folder = Path(reconstituted_path, part)
        if not next_folder.exists():
            error = (
                reconstituted_path,
                fuzzy_match([part], os.listdir(reconstituted_path))[0][1],
                part,
            )
            break
        reconstituted_path = str(next_folder)

    msg = (
        f"The following file could not be found: '{path}'. "
        f"In the folder '{error[0]}' the following '{error[1]}' "
        f"is the closest match to the supplied '{error[2]}', did "
        "you mean to open this?"
    )
    raise ValueError(msg)


def fuzzy_match(
    user_input: Iterable[str], possibilities: Iterable[str]
) -> list[tuple[str, str]]:
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
    results: list[tuple[str, str]] = list()

    for user_word in transform_to_list(user_input):
        best: tuple[float, Any] = (-1, "")
        for possibility in possibilities:
            distance = difflib.SequenceMatcher(
                a=user_word, b=possibility
            ).ratio()  # noqa: E501
            if distance > best[0]:
                best = (distance, possibility)
        results.append((user_word, best[1]))

    return results


def update_step_contents_to_step_names(
    prev_names: Iterable[str], new_names: Iterable[str], folder: FilePath
) -> None:
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


def update_step_names_in_subfolders(
    folder: Path, prev_names: Iterable[str], new_names: Iterable[str]
) -> None:
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


def update_step_names_in_file(
    file_: Path, prev_names: Iterable[str], new_names: Iterable[str]
) -> None:
    """Update step names in file following the `--restart` option."""
    try:
        text = file_.read_text()
    except UnicodeDecodeError as err:
        log.warning(f"Failed to read file {file_}. Error is {err}")
        return
    for s1, s2 in zip(prev_names, new_names):
        text = text.replace(s1, s2)
    file_.write_text(text)
    return


def check_CNS_usage(modules_params: ParamMap) -> None:
    """Check that a topology module is run prior to modules requiring CNS.

    Parameters
    ----------
    modules_params : ParamMap
        Dict of modules parameters.
        Only used to obtain ordered list of modules.

    Raises
    ------
    DependencyError
        Error thrown if topology not run before a CNS module.
    """
    from haddock.core.defaults import CNS_MODULES
    generated_topology: bool = False
    for _module_name in modules_params:
        module_name = get_module_name(_module_name)
        # Check if this module is a topology module
        if modules_category[module_name] == "topology":
            # Set the flag
            generated_topology = True
        # Check if this module is a CNS module (that require topology)
        if module_name in CNS_MODULES:
            # Check that topology was generated
            if not generated_topology:
                raise DependencyError(
                    msg="A topology module should be used prior to CNS module.",
                    module=module_name,
                    dependency=", ".join([
                        k for k, v in modules_category.items()
                        if v == "topology"
                        ])
                    )
            # We can stop here as either error raised or check passsed
            break
