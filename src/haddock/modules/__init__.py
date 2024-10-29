"""HADDOCK3 modules."""

import re

from abc import ABC, abstractmethod
from contextlib import contextmanager, suppress
from copy import deepcopy
from functools import partial
from os import linesep
from pathlib import Path

from haddock import EmptyPath, log, modules_defaults_path
from haddock.core.defaults import MODULE_IO_FILE, INTERACTIVE_RE_SUFFIX
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import (
    Any,
    Container,
    FilePath,
    Generator,
    Literal,
    Optional,
    ParamDict,
    Union,
)
from haddock.gear import config
from haddock.gear.clean_steps import clean_output
from haddock.gear.known_cns_errors import find_all_cns_errors
from haddock.gear.parameters import config_mandatory_general_parameters
from haddock.gear.yaml2cfg import read_from_yaml_config, find_incompatible_parameters
from haddock.libs.libhpc import HPCScheduler
from haddock.libs.libio import folder_exists, working_directory
from haddock.libs.libmpi import MPIScheduler
from haddock.libs.libontology import ModuleIO, PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libtimer import log_time
from haddock.libs.libutil import recursive_dict_update


modules_folder = Path(__file__).resolve().parent

_folder_match_regex = "[a-zA-Z]*/"
modules_category = {
    module.name: category.name
    for category in modules_folder.glob(_folder_match_regex)
    for module in category.glob(_folder_match_regex)
}
"""Indexes each module in its specific category. Keys are Paths to the module,
values are their categories. Categories are the modules parent folders."""

modules_names = set(modules_category.keys())

category_hierarchy = [
    "topology",
    "sampling",
    "refinement",
    "scoring",
    "analysis",
    "extras",
]

# this dictionary defines non-mandatory general parameters that can be defined
# as global parameters thus affect all modules, or, instead, can be defined per
# module where the module definition overwrites global definition. Not all
# modules will use these parameters. It is the responsibility of the module to
# extract the parameters it needs.
# the config file is in modules/defaults.cfg
non_mandatory_general_parameters_defaults = read_from_yaml_config(
    modules_defaults_path
)  # noqa : E501

incompatible_defaults_params = find_incompatible_parameters(modules_defaults_path)

config_readers = {
    ".yaml": read_from_yaml_config,
    ".cfg": config.load,
}

_step_folder_regex = tuple(
    r"[0-9]+_" + mod_name for mod_name in modules_category.keys()
)
step_folder_regex = "(" + "|".join(_step_folder_regex) + ")"
"""
String for regular expression to match module folders in a run directory.

It will match folders with a numeric prefix followed by underscore ("_")
followed by the name of a module.

Example: https://regex101.com/r/roHls9/1
"""

step_folder_regex_re = re.compile(step_folder_regex)
"""
Compiled regular expression from :py:const:`step_folder_regex`.

It will match folders with a numeric prefix followed by underscore ("_")
followed by the name of a module.

Example: https://regex101.com/r/roHls9/1
"""


@contextmanager
def _not_valid_config() -> Generator[None, None, None]:
    try:
        yield
    except KeyError as err:
        emsg = (
            "The configuration file extension is not supported. "
            f"Supported types are {', '.join(config_readers.keys())}."
        )
        raise ConfigurationError(emsg) from err


class BaseHaddockModule(ABC):
    """HADDOCK3 module's base class."""

    name: str

    def __init__(self, order: int, path: Path, params_fname: FilePath) -> None:
        """
        HADDOCK3 modules base class.

        Parameters
        ----------
        params : dict or path to HADDOCK3 configuration file
            A dictionary or a path to a HADDOCK3 configuration file
            containing the initial module parameters. Usually this is
            defined by the default params.
        """
        self.order = order
        self.path = path
        self.previous_io = self._load_previous_io()

        # instantiate module's parameters
        self._origignal_config_file = params_fname
        with _not_valid_config():
            extension = Path(params_fname).suffix
            self._original_params = config_readers[extension](params_fname)

        self._params: ParamDict = {}
        self.update_params(update_from_cfg_file=params_fname)

    @property
    def params(self) -> ParamDict:
        """Configuration parameters."""  # noqa: D401
        return self._params

    def reset_params(self) -> None:
        """Reset parameters to the ones used to instantiate the class."""
        self._params.clear()
        self.update_params(**self._original_params)

    def update_params(
        self,
        update_from_cfg_file: Optional[FilePath] = None,
        **params: Any,
    ) -> None:
        """
        Update the modules parameters.

        Add/update to the current modules parameters the ones given in
        the function call. If you want to enterily replace the modules
        parameters to their default values use the `reset_params()`
        method.

        Update takes places recursively, that is, nested dictionaries
        will be updated accordingly.

        To update the current config with the parameters defined in an
        HADDOCK3 configuration file use the `update_from_cfg_file`
        parameter.

        To update from a JSON file, first load the JSON into a
        dictionary and unpack the dictionary to the function call.

        Examples
        --------
        >>> m.update_params(param1=value1, param2=value2)

        >>> m.update_params(**param_dict)

        >>> m.update_params(update_from_cfg_file=path_to_file)

        # if you wish to start from scratch
        >>> m.reset_params()
        >>> m.update_params(...)
        """
        if update_from_cfg_file and params:
            _msg = (
                "You can not provide both `update_from_cfg_file` " "and key arguments."
            )
            raise TypeError(_msg)

        if update_from_cfg_file:
            with _not_valid_config():
                extension = Path(update_from_cfg_file).suffix
                params = config_readers[extension](update_from_cfg_file)

        # the updating order is relevant
        _n = recursive_dict_update(
            non_mandatory_general_parameters_defaults, self._params
        )
        self._params = recursive_dict_update(_n, params)
        self._fill_emptypaths()
        self._confirm_fnames_exist()

    def save_config(self, path: FilePath) -> None:
        """Save current parameters to a HADDOCK3 config file."""
        # creates this dictionary for the config to have the module name
        # key in brackets, for example:
        #
        # [topoaa]
        # ...
        ignore = config_mandatory_general_parameters.union(
            non_mandatory_general_parameters_defaults
        )  # noqa: 501
        params = deepcopy(self.params)

        with suppress(KeyError):
            for key in list(ignore):
                params.pop(key)

        config.save({self.name: params}, path)

    def add_parent_to_paths(self) -> None:
        """Add parent path to paths."""
        # convert paths to relative by appending parent
        for key, value in self.params.items():
            if value and key.endswith("_fname"):
                if not Path(value).is_absolute():
                    self.params[key] = Path("..", value)
        return

    @abstractmethod
    def _run(self) -> None: ...

    def run(self, **params: Any) -> None:
        """Execute the module."""
        log.info(f"Running [{self.name}] module")

        self.update_params(**params)
        self.add_parent_to_paths()

        with working_directory(self.path):
            self._run()

        log.info(f"Module [{self.name}] finished.")

    def clean_output(self) -> None:
        """
        Clean module output folder.

        See Also
        --------
        :py:func:`haddock.gear.clean_steps.clean_output`
        """
        with log_time("cleaning output files took"):
            clean_output(self.path, self.params["ncores"])

    @classmethod
    @abstractmethod
    def confirm_installation(cls) -> None:
        """
        Confirm the third-party software needed for the module is installed.

        HADDOCK3's own modules should just return.
        """
        return

    def export_io_models(self, faulty_tolerance: float = 0.0) -> None:
        """
        Export input/output to the ModuleIO interface.

        Modules that do not perform any operation on PDB files should have
         input = output.

        This function implements a common interface for all modules.

        Parameters
        ----------
        faulty_tolerance : int, default 0
            The percentage of missing output allowed. If 20 is given,
            raises an error if 20% of the expected output is missing (not
            saved to disk).
        """
        self.output_models: Union[list[PDBFile], dict[int, PDBFile]]
        assert self.output_models, "`self.output_models` cannot be empty."
        io = ModuleIO()
        # add the input models
        io.add(self.previous_io.output, "i")
        # add the output models
        io.add(self.output_models, "o")
        # Removes un-generated outputs and compute percentage of ungenerated
        faulty = io.check_faulty()
        # Save outputs
        io.save()
        # Check if number of generated outputs is under the tolerance threshold
        if faulty > faulty_tolerance:
            _msg = (
                f"{faulty:.2f}% of output was not generated for this module "
                f"and tolerance was set to {faulty_tolerance:.2f}%."
            )
            # Try to detect CNS errors
            if detected_errors := find_all_cns_errors(self.path):
                _msg += linesep
                for error in detected_errors.values():
                    _msg += f'{str(error["error"])}{linesep}'
            # Show final error message
            self.finish_with_error(_msg)

    def finish_with_error(self, reason: object = "Module has failed.") -> None:
        """Finish with error message."""
        if isinstance(reason, Exception):
            raise RuntimeError("Module has failed.") from reason

        else:
            raise RuntimeError(reason)

    def _load_previous_io(
        self,
        filename: FilePath = MODULE_IO_FILE,
    ) -> ModuleIO:
        if self.order == 0:
            self._num_of_input_molecules = 0
            return ModuleIO()

        io = ModuleIO()
        previous_io = Path(self.previous_path(), filename)

        if previous_io.is_file():
            io.load(previous_io)

        self._num_of_input_molecules = len(io.output)

        return io

    def previous_path(self) -> Path:
        """Give the path from the previous calculation."""
        previous = get_module_steps_folders(self.path.resolve().parent)

        try:
            # return Path(previous[self.order - 1])
            return self.last_step_folder(previous, self.order - 1)
        except IndexError:
            return self.path

    @staticmethod
    def last_step_folder(folders, index):
        """Retrieve last step folder."""
        with_ind = [folder for folder in folders if int(folder.split("_")[0]) == index]
        nb_with_ind = len(with_ind)
        # No matching index
        if nb_with_ind == 0:
            raise IndexError
        # Only one matching index
        elif nb_with_ind == 1:
            return with_ind[0]
        # Case of multiple matching index
        else:
            for folder in with_ind:
                if folder.split("_")[-1] != INTERACTIVE_RE_SUFFIX:
                    return folder
            return with_ind[0]

    def log(self, msg: str, level: str = "info") -> None:
        """
        Log a message with a common header.

        Currently the header is the [MODULE NAME] in square brackets.

        Parameters
        ----------
        msg : str
            The log message.

        level : str
            The level log: 'debug', 'info', ...
            Defaults to 'info'.
        """
        getattr(log, level)(f"[{self.name}] {msg}")

    def _confirm_fnames_exist(self) -> None:
        for param, value in self._params.items():
            if param.endswith("_fname") and value:
                if not Path(value).exists():
                    raise FileNotFoundError(f"File not found: {str(value)!r}")

    def _fill_emptypaths(self) -> None:
        """Fill empty paths."""
        for param, value in list(self._params.items()):
            if param.endswith("_fname") and not value:
                self._params[param] = EmptyPath()


EngineMode = Literal["batch", "local", "mpi"]


def get_engine(
    mode: str,
    params: dict[Any, Any],
) -> partial[Union[HPCScheduler, Scheduler, MPIScheduler]]:
    """
    Create an engine to run the jobs.

    Parameters
    ----------
    mode : str
        The type of engine to create

    params : dict
        A dictionary containing parameters for the engine.
        `get_engine` will retrieve from `params` only those parameters
        needed and ignore the others.
    """
    # a bit of a factory pattern here
    # this might end up in another module but for now its fine here
    if mode == "batch":
        return partial(  # type: ignore
            HPCScheduler,
            target_queue=params["queue"],
            queue_limit=params["queue_limit"],
            concat=params["concat"],
        )

    elif mode == "local":
        return partial(  # type: ignore
            Scheduler,
            ncores=params["ncores"],
            max_cpus=params["max_cpus"],
        )
    elif mode == "mpi":
        return partial(MPIScheduler, ncores=params["ncores"])  # type: ignore

    else:
        available_engines = ("batch", "local", "mpi")
        raise ValueError(
            f"Scheduler `mode` {mode!r} not recognized. "
            f"Available options are {', '.join(available_engines)}"
        )


def get_module_steps_folders(
    folder: FilePath,
    modules: Optional[Container[int]] = None,
) -> list[str]:
    """
    Return a sorted list of the step folders in a running directory.

    Example
    -------
    Consider the folder structure:

    run_dir/
        0_topoaa/
        1_rigidbody/
        2_caprieval/
        3_bad_module_name/
        data/

    >>> get_module_steps_folders("run_dir")
    >>> ["0_topoaa", "1_rigidbody", "2_caprieval"]

    Parameters
    ----------
    folder : str or Path
        Path to the run directory, or to the folder containing the step
        folders.

    Returns
    -------
    list of str
        List containing strings with the names of the step folders.
    """
    folders = (p.name for p in Path(folder).iterdir() if p.is_dir())
    steps = sorted(
        (f for f in folders if step_folder_regex_re.search(f)),
        key=lambda x: int(x.split("_")[0]),
    )
    if modules:
        steps = [
            st
            for st in steps
            if all(
                [
                    int(st.split("_")[0]) in modules,
                    st.split("_")[1] in modules_names,
                ]
            )
        ]
    return steps


def is_step_folder(path: FilePath) -> bool:
    """
    Assess whether a folder is a possible step folder.

    The folder is considered a step folder if has a zero or positive
    integer index followed by a name of a module.

    Parameters
    ----------
    path : str or pathlib.Path
        The path to the folder.

    Returns
    -------
    bool
        Whether the folder is a step folder or not.
    """
    path = Path(path)
    folder_exists(path)
    main_folder_name = path.name
    parts = main_folder_name.split("_")
    if len(parts) == 2 and parts[0].isdigit() and parts[1] in modules_category:
        return True
    else:
        return False
