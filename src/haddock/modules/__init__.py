"""Workflow module logic"""
import os
import logging
import contextlib
from pathlib import Path
import toml
from haddock.error import StepError
from haddock.ontology import ModuleIO
from haddock.defaults import MODULE_PATH_NAME, MODULE_IO_FILE, TOPOLOGY_PATH


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


class BaseHaddockModule:
    """Base class for any HADDOCK module"""
    def __init__(self, order, path, cns_script="", defaults=""):
        self.order = order
        self.path = path
        self.previous_io = self._load_previous_io()

        if cns_script:
            self.cns_folder_path = cns_script.resolve().parent
            self.cns_protocol_path = cns_script

        self.defaults = defaults

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
    def defaults(self):
        return self._defaults

    @defaults.setter
    def defaults(self, path):
        try:
            self._defaults = toml.load(path)
        except FileNotFoundError as err:
            _msg = f"Default configuration file path not found: {str(path)!r}"
            raise StepError(_msg) from err

    def run(self, params):
        self.update_defaults(**params)

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
            return previous[-2]
        except IndexError:
            return self.path

    def update_defaults(self, **parameters):
        """Update defaults parameters with run-specific parameters."""
        self._defaults.update(parameters)


@contextlib.contextmanager
def working_directory(path):
    """Changes working directory and returns to previous on exit"""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
