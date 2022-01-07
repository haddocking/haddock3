"""HADDOCK3 modules."""
from abc import ABC, abstractmethod
from copy import deepcopy
from functools import partial
from pathlib import Path

from haddock import log as log
from haddock.core.defaults import MODULE_IO_FILE, cns_exec
from haddock.gear.config_reader import read_config
from haddock.libs.libhpc import (
    HPCScheduler,
    HPCScheduler_CONCAT_DEFAULT,
    HPCWorker_QUEUE_DEFAULT,
    HPCWorker_QUEUE_LIMIT_DEFAULT,
    )
from haddock.libs.libio import working_directory
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


# non-mandatory general parameters that can be defined as global parameters
# that affect all modules, or given per module where local overwrites global.
# Not all modules will use these parameters. It is the responsibility of the
# module to extract those parameters it needs.
non_mandatory_general_parameters_defaults = {
    "concat": HPCScheduler_CONCAT_DEFAULT,
    "cns_exec": cns_exec,
    "mode": "local",
    "ncores": 8,
    "queue": HPCWorker_QUEUE_DEFAULT,
    "queue_limit": HPCWorker_QUEUE_LIMIT_DEFAULT,
    "self_contained": False,
    }


class BaseHaddockModule(ABC):
    """HADDOCK3 module's base class."""

    def __init__(self, order, path, params):
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
        self.params = params

    @property
    def params(self):
        """Configuration parameters."""  # noqa: D401
        return self._params

    @params.setter
    def params(self, path_or_dict):
        if isinstance(path_or_dict, Path):
            conf_dict = read_config(path_or_dict)
        elif isinstance(path_or_dict, dict):
            conf_dict = path_or_dict
        else:
            _msg = (
                "Argument does not satisfy condition, must be path or "
                f"dict. {type(path_or_dict)} given."
                )
            raise TypeError(_msg)

        # the new parameters are created on top of the default general
        # parameters for modules because the general parameters for modules are
        # not module specific, instead they are used by several modules.
        _d = deepcopy(non_mandatory_general_parameters_defaults)
        self._params = recursive_dict_update(_d, conf_dict)

        for param, value in self._params.items():
            if param.endswith('_fname'):
                if not Path(value).exists():
                    raise FileNotFoundError(f'File not found: {str(value)!r}')

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

        self.params = params
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

    def finish_with_error(self, reason="Module has failed."):
        """Finish with error message."""
        if isinstance(reason, Exception):
            raise RuntimeError("Module has failed.") from reason

        else:
            raise RuntimeError(reason)

    def _load_previous_io(self):
        if self.order == 0:
            return ModuleIO()

        io = ModuleIO()
        previous_io = self.previous_path() / MODULE_IO_FILE
        if previous_io.is_file():
            io.load(previous_io)
        return io

    def previous_path(self):
        """Give the path from the previous calculation."""
        previous = sorted(list(self.path.resolve().parent.glob('[0-9][0-9]*/')))
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

    else:
        available_engines = ('hpc', 'local')
        raise ValueError(
            f"Scheduler `mode` {mode!r} not recognized. "
            f"Available options are {', '.join(available_engines)}"
            )
