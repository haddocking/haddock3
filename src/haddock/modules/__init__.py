"""HADDOCK3 modules."""
import re
from abc import ABC, abstractmethod
from contextlib import contextmanager
from functools import partial
from pathlib import Path

from haddock import EmptyPath, log, modules_defaults_path
from haddock.core.defaults import MODULE_IO_FILE
from haddock.core.exceptions import ConfigurationError
from haddock.gear.config_reader import read_config
from haddock.gear.config_writer import convert_config as _convert_config
from haddock.gear.config_writer import save_config as _save_config
from haddock.gear.parameters import config_mandatory_general_parameters
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs.libhpc import HPCScheduler
from haddock.libs.libio import working_directory
from haddock.libs.libmpi import MPIScheduler
from haddock.libs.libontology import ModuleIO
from haddock.libs.libparallel import Scheduler
from haddock.libs.libutil import recursive_dict_update


modules_folder = Path(__file__).resolve().parent

_folder_match_regex = '[a-zA-Z]*/'
modules_category = {
    module.name: category.name
    for category in modules_folder.glob(_folder_match_regex)
    for module in category.glob(_folder_match_regex)
    }
"""Indexes each module in its specific category. Keys are Paths to the module,
values are their categories. Categories are the modules parent folders."""

category_hierarchy = [
    "topology",
    "sampling",
    "refinement",
    "scoring",
    "analysis",
    ]

# this dictionary defines non-mandatory general parameters that can be defined
# as global parameters thus affect all modules, or, instead, can be defined per
# module where the module definition overwrites global definition. Not all
# modules will use these parameters. It is the responsibility of the module to
# extract the parameters it needs.
# the config file is in modules/defaults.cfg
non_mandatory_general_parameters_defaults = \
    read_from_yaml_config(modules_defaults_path)

config_readers = {
    ".yaml": read_from_yaml_config,
    ".cfg": read_config,
    }

_step_folder_regex = tuple(
    r"[0-9]+_" + mod_name
    for mod_name in modules_category.keys()
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
def _not_valid_config():
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

    def __init__(self, order, path, params_fname):
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

        self._params = {}
        self.update_params(update_from_cfg_file=params_fname)

    @property
    def params(self):
        """Configuration parameters."""  # noqa: D401
        return self._params

    def reset_params(self):
        """Reset parameters to the ones used to instantiate the class."""
        self._params.clear()
        self.update_params(**self._original_params)

    def update_params(self, update_from_cfg_file=None, **params):
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
                "You can not provide both `update_from_cfg_file` "
                "and key arguments."
                )
            raise TypeError(_msg)

        if update_from_cfg_file:
            with _not_valid_config():
                extension = Path(update_from_cfg_file).suffix
                params = config_readers[extension](update_from_cfg_file)

        # the updating order is relevant
        _n = recursive_dict_update(
            non_mandatory_general_parameters_defaults,
            self._params)
        self._params = recursive_dict_update(_n, params)
        self._fill_emptypaths()
        self._confirm_fnames_exist()

    def save_config(self, path):
        """Save current parameters to a HADDOCK3 config file."""
        # creates this dictionary for the config to have the module name
        # key in brackets, for example:
        #
        # [topoaa]
        # ...
        save_config_ignored({self.name: self.params}, path)

    def add_parent_to_paths(self):
        """Add parent path to paths."""
        # convert paths to relative by appending parent
        for key, value in self.params.items():
            if value and key.endswith('_fname'):
                if not Path(value).is_absolute():
                    self.params[key] = Path('..', value)
        return

    def run(self, **params):
        """Execute the module."""
        log.info(f'Running [{self.name}] module')

        self.update_params(**params)
        self.add_parent_to_paths()

        with working_directory(self.path):
            self._run()

        log.info(f'Module [{self.name}] finished.')

    @classmethod
    @abstractmethod
    def confirm_installation(self):
        """
        Confirm the third-party software needed for the module is installed.

        HADDOCK3's own modules should just return.
        """
        return

    def export_output_models(self, faulty_tolerance=0):
        """
        Export output to the ModuleIO interface.

        Modules that generate PDBs that other models should take as input,
        should export those PDBs registries through the ModuleIO interface.

        This function implements a common interface for all modules requiring
        this feature.

        Parameters
        ----------
        faulty_tolerance : int, default 0
            The percentage of missing output allowed. If 20 is given,
            raises an error if 20% of the expected output is missing (not
            saved to disk).
        """
        assert self.output_models, "`self.output_models` cannot be empty."
        io = ModuleIO()
        io.add(self.output_models, "o")
        faulty = io.check_faulty()
        if faulty > faulty_tolerance:
            _msg = (
                f"{faulty:.2f}% of output was not generated for this module "
                f"and tolerance was set to {faulty_tolerance:.2f}%.")
            self.finish_with_error(_msg)
        io.save()

    def finish_with_error(self, reason="Module has failed."):
        """Finish with error message."""
        if isinstance(reason, Exception):
            raise RuntimeError("Module has failed.") from reason

        else:
            raise RuntimeError(reason)

    def _load_previous_io(self, filename=MODULE_IO_FILE):
        if self.order == 0:
            self._num_of_input_molecules = 0
            return ModuleIO()

        io = ModuleIO()
        previous_io = Path(self.previous_path(), filename)

        if previous_io.is_file():
            io.load(previous_io)

        self._num_of_input_molecules = len(io.output)

        return io

    def previous_path(self):
        """Give the path from the previous calculation."""
        previous = get_module_steps_folders(self.path.resolve().parent)

        try:
            return Path(previous[self.order - 1])
        except IndexError:
            return self.path

    def log(self, msg, level='info'):
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
        getattr(log, level)(f'[{self.name}] {msg}')

    def _confirm_fnames_exist(self):
        for param, value in self._params.items():
            if param.endswith('_fname') and value:
                if not Path(value).exists():
                    raise FileNotFoundError(f'File not found: {str(value)!r}')

    def _fill_emptypaths(self):
        """Fill empty paths."""
        for param, value in list(self._params.items()):
            if param.endswith('_fname') and not value:
                self._params[param] = EmptyPath()


def get_engine(mode, params):
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
    if mode == 'hpc':
        return partial(
            HPCScheduler,
            target_queue=params['queue'],
            queue_limit=params['queue_limit'],
            concat=params['concat'],
            )

    elif mode == 'local':
        return partial(
            Scheduler,
            ncores=params['ncores'],
            )
    elif mode == "mpi":
        return partial(MPIScheduler, ncores=params["ncores"])

    else:
        available_engines = ("hpc", "local", "mpi")
        raise ValueError(
            f"Scheduler `mode` {mode!r} not recognized. "
            f"Available options are {', '.join(available_engines)}"
            )


def convert_config(params):
    """
    Convert a module's parameters dictionary to a HADDOCK3 user config text.

    This function is a generator.

    Examples
    --------
    >>> gen = convert_config(params)
    >>> text = os.linesep.join(gen)
    >>> with open("params.cfg", "w") as fout:
    >>>     fout.write(text)

    Parameters
    ----------
    params : dictionary
        The dictionary containing the parameters.

    Yields
    ------
    str
        Line by line for the HADDOCK3 user configuration file.

    See Also
    --------
    :py:func:`haddock.gear.config_writer.convert_config`.
    """
    return _convert_config(
        params,
        ignore_params=non_mandatory_general_parameters_defaults,
        module_names=set(modules_category.keys()),
        )


def save_config(*args, **kwargs):
    """
    Save HADDOCK3 configuration dictionary to user config file.

    Saves all parameters, even the `non_mandatory_general_parameters_defaults`
    which can be repeated between the module parameters and the general
    parameters.

    Parameters
    ----------
    params : dict
        The dictionary containing the parameters.

    path : str or pathlib.Path
        File name where to save the configuration file.

    See Also
    --------
    :py:func:`haddock.gear.config_writer.save_config`.
    """
    kwargs.setdefault("module_names", set(modules_category.keys()))
    return _save_config(*args, **kwargs)


def save_config_ignored(*args, **kwargs):
    """
    Save HADDOCK3 configuration dictionary to user config file.

    Ignores the
    :py:data:`haddock.modules.non_mandatory_general_parameters_defaults`
    parameters.

    Useful to keep clean versions of the modules' specific parameters.

    Parameters
    ----------
    params : dict
        The dictionary containing the parameters.

    path : str or pathlib.Path
        File name where to save the configuration file.

    See Also
    --------
    :py:func:`save_config`
    :py:func:`haddock.gear.config_writer.save_config`.
    """
    kwargs.setdefault("module_names", set(modules_category.keys()))
    kwargs.setdefault(
        "ignore_params",
        config_mandatory_general_parameters.union(non_mandatory_general_parameters_defaults),  # noqa: 501
        )
    return _save_config(*args, **kwargs)


def get_module_steps_folders(folder):
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
    return steps
