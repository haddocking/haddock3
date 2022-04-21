"""HADDOCK3 modules."""
from abc import ABC, abstractmethod
from contextlib import contextmanager
from functools import partial
from pathlib import Path

from haddock import EmptyPath, log, modules_defaults_path
from haddock.core.defaults import MODULE_IO_FILE
from haddock.core.exceptions import ConfigurationError
from haddock.gear.config_reader import read_config
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
        self._params = {}
        self.update_params(update_from_cfg_file=params_fname)

    @property
    def params(self):
        """Configuration parameters."""  # noqa: D401
        return self._params

    def reset_params(self):
        """Reset parameters to the ones used to instantiate the class."""
        self._params.clear()
        self.update_params(update_from_cfg_file=self._original_config_file)

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
            return ModuleIO()
        io = ModuleIO()
        previous_io = self.previous_path() / filename
        if previous_io.is_file():
            io.load(previous_io)
        return io

    def previous_path(self):
        """Give the path from the previous calculation."""
        # [0-9]* below is not a regex, is a bash wildkey
        previous = sorted(list(self.path.resolve().parent.glob('[0-9]*/')))
        try:
            return previous[self.order - 1]
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
