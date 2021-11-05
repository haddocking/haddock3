"""Workflow module logic"""
import os
import logging
import contextlib
from abc import ABC, abstractmethod
from pathlib import Path

from haddock.core.defaults import MODULE_PATH_NAME, MODULE_IO_FILE
from haddock.core.exceptions import StepError
from haddock.gear.config_reader import read_config
from haddock.libs.libontology import ModuleIO


logger = logging.getLogger(__name__)

modules_folder = Path(__file__).resolve().parent

_folder_match_regex = '[a-zA-Z]*/'
modules_category = {
    module.name: category.name
    for category in modules_folder.glob(_folder_match_regex)
    for module in category.glob(_folder_match_regex)
    }
"""Indexes each module in its specific category. Keys are Paths to the module,
values are their categories. Categories are the modules parent folders."""


general_parameters_affecting_modules = {'ncores', 'cns_exec'}
"""These parameters are general parameters that may be applicable to modules
specifically. Therefore, they should be considered as part of the "default"
module's parameters. Usually, this set is used to filter parameters during
the run prepraration phase. See, `gear.prepare_run`."""


class BaseHaddockModule(ABC):
    def __init__(self, order, path, params, cns_script=""):
        """
        Base class for any HADDOCK module

        Parameters
        ----------
        params : dict or path to HADDOCK3 configuration file
            A dictionary or a path to an HADDOCK3 configuration file
            containing the initial module parameters. Usually this is
            defined by the default params.
        """
        self.order = order
        self.path = path
        self.previous_io = self._load_previous_io()

        if cns_script:
            self.cns_folder_path = cns_script.resolve().parent
            self.cns_protocol_path = cns_script

        self.params = params

        try:
            with open(self.cns_protocol_path) as input_handler:
                self.recipe_str = input_handler.read()
        except FileNotFoundError:
            _msg = f"Error while opening workflow {self.cns_protocol_path}"
            raise StepError(_msg)
        except AttributeError:
            # No CNS-like module
            pass

    @property
    def params(self):
        return self._params

    @params.setter
    def params(self, path_or_dict):
        if isinstance(path_or_dict, dict):
            self._params = path_or_dict
        else:
            try:
                self._params = read_config(path_or_dict)
            except FileNotFoundError as err:
                _msg = f"Default configuration file not found: {str(path_or_dict)!r}"
                raise FileNotFoundError(_msg) from err
            except TypeError as err:
                _msg = (
                    "Argument does not satisfy condition, must be path or "
                    f"dict. {type(path_or_dict)} given."
                    )
                raise TypeError(_msg) from err

    @abstractmethod
    def run(self, params):
        self.update_params(**params)
        self.params.setdefault('ncores', None)
        self.params.setdefault('cns_exec', None)

    @classmethod
    @abstractmethod
    def confirm_installation(self):
        """
        Confirm the third-party software needed for the module is installed.

        HADDOCK3's own modules should just return.
        """
        return

    def finish_with_error(self, message=""):
        if not message:
            message = "Module has failed"
        logger.error(message)
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
        previous = sorted(list(self.path.resolve().parent.glob('[0-9][0-9]*/')))
        try:
            return previous[self.order - 1]
        except IndexError:
            return self.path

    def update_params(self, **parameters):
        """Update defaults parameters with run-specific parameters."""
        self._params.update(parameters)


@contextlib.contextmanager
def working_directory(path):
    """Changes working directory and returns to previous on exit"""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
