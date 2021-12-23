"""Functionalities related to CNS modules."""
import os
import shutil
from pathlib import Path

from haddock import toppar_path as global_toppar
from haddock.modules import BaseHaddockModule


class BaseCNSModule(BaseHaddockModule):
    """
    Operation module for CNS.

    Contains additional functionalities excusive for CNS modules.
    """

    def __init__(self, *args, cns_script="", **kwargs):
        """
        Instantitate CNSModule.

        Parameters
        ----------
        cns_script : str or pathlib.Path
            Path to the main module's cns script.
        """
        super().__init__(*args, **kwargs)

        self.cns_folder_path = Path(cns_script.resolve().parent)
        self.cns_protocol_path = cns_script
        self.toppar_path = global_toppar
        self.recipe_str = self.cns_protocol_path.read_text()

    def run(self, **params):
        """Execute the module."""
        self.envvars = self.default_envvars()

        if params['self_contained']:
            self.make_self_contained()

        super().run(**params)

    def default_envvars(self):
        """Return default env vars updated to `envvars` (if given)."""
        default_envvars = {
            "MODULE": str(self.cns_folder_path),
            "MODDIR": ".",
            "TOPPAR": str(self.toppar_path),
            }
        print(default_envvars)

        return default_envvars

    def save_envvars(self, filename="envvars"):
        """Save envvars needed for CNS to a file in the module's folder."""
        # there are so few variables, best to handle them by hand
        lines = (
            "#!/bin/bash",
            "# for debugging purposes source this file from within the ",
            "# module folder for example, from within '00_topoaa'",
            "export MODULE=cns",
            "export MODDIR=.",
            "export TOPPAR=../toppar",
            )

        fstr = os.linesep.join(lines)
        Path(self.path, filename).write_text(fstr)
        return

    def make_self_contained(self):
        """Create folders to make run self-contained."""
        _ = Path(self.path, "cns")
        shutil.copytree(self.cns_folder_path, _)
        self.cns_folder_path = Path(".", "cns")

        self.cns_protocol_path = Path(
            self.cns_folder_path,
            self.cns_protocol_path.name,
            )

        if not Path(self.toppar_path.name).exists():
            shutil.copytree(self.toppar_path, self.toppar_path.name)
        self.toppar_path = Path("..", self.toppar_path.name)

        self.envvars = self.default_envvars()
        self.save_envvars()
