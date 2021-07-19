"""Workflow module logic"""
import os
import logging
import contextlib
from pathlib import Path
import toml
from haddock.error import RecipeError
from haddock.ontology import ModuleIO
from haddock.defaults import MODULE_PATH_NAME, MODULE_IO_FILE, TOPOLOGY_PATH

logger = logging.getLogger(__name__)


class BaseHaddockModule:
    """Base class for any HADDOCK module"""
    def __init__(self, order, path, cns_script="", defaults=""):
        self.order = order
        self.path = path
        self.previous_io = self._load_previous_io()

        if cns_script:
            self.cns_folder_path = cns_script.resolve().parent.absolute()
            self.cns_recipe_path = cns_script
        if defaults:
            self.defaults_path = defaults

        try:
            with open(self.cns_recipe_path) as input_handler:
                self.recipe_str = input_handler.read()
        except FileNotFoundError:
            raise RecipeError(f"Error while opening recipe {self.cns_recipe_path}")
        except AttributeError:
            # No CNS-like module
            pass
        try:
            self.defaults = toml.load(self.defaults_path)
        except FileNotFoundError:
            raise RecipeError(f"Error while opening defaults {self.defaults_path}")
        except AttributeError:
            # No CNS-like module
            pass

    def run(self, module_information):
        raise NotImplementedError()

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
        if self.order > 1:
            return self.path.resolve().parent.absolute() / f"{MODULE_PATH_NAME}{self.order-1}"
        if self.order == 1:
            return self.path.resolve().parent.absolute() / TOPOLOGY_PATH

        return self.path

    def patch_defaults(self, module_parameters):
        """Apply custom module parameters given to defaults dictionary"""
        for k in module_parameters:
            if k in self.defaults["params"]:
                self.defaults["params"][k] = module_parameters[k]


@contextlib.contextmanager
def working_directory(path):
    """Changes working directory and returns to previous on exit"""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
