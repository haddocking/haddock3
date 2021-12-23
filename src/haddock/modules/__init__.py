"""HADDOCK3 modules."""
from abc import ABC, abstractmethod
from functools import partial
from pathlib import Path

from haddock import log as log
from haddock.core.defaults import MODULE_IO_FILE
from haddock.gear.config_reader import read_config
from haddock.libs.libio import working_directory
from haddock.libs.libhpc import HPCScheduler
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


general_parameters_affecting_modules = {
    'cns_exec',
    'concat',
    'mode',
    'ncores',
    'queue',
    'queue_limit',
    'self_contained',
    }
"""These parameters are general parameters that may be applicable to modules
specifically. Therefore, they should be considered as part of the "default"
module's parameters. Usually, this set is used to filter parameters during
the run prepraration phase. See, `gear.prepare_run`."""


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
        if isinstance(path_or_dict, dict):
            self._params = path_or_dict
        else:
            try:
                self._params = read_config(path_or_dict)
            except FileNotFoundError as err:
                _msg = (
                    "Default configuration file not found: "
                    f"{str(path_or_dict)!r}"
                    )
                raise FileNotFoundError(_msg) from err
            except TypeError as err:
                _msg = (
                    "Argument does not satisfy condition, must be path or "
                    f"dict. {type(path_or_dict)} given."
                    )
                raise TypeError(_msg) from err

        for param, value in self._params.items():
            if param.endswith('_fname'):
                if not Path(value).exists():
                    raise FileNotFoundError(f'File not found: {str(value)!r}')

    def run(self, **params):
        """Execute the module."""
        log.info(f'Running [{self.name}] module')
        self.update_params(**params)
        self.params.setdefault('ncores', None)
        self.params.setdefault('cns_exec', None)
        self.params.setdefault('mode', None)
        self.params.setdefault('concat', None)
        self.params.setdefault('queue_limit', None)

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

    def finish_with_error(self, message=""):
        """Finish with error message."""
        if not message:
            message = "Module has failed"
        log.error(message)
        raise SystemExit

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

    def update_params(self, **parameters):
        """Update defaults parameters with run-specific parameters."""
        self.params = recursive_dict_update(self._params, parameters)

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
