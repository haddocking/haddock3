"""Functionalities related to CNS modules."""
import os
import shutil
from pathlib import Path

from haddock import log
from haddock import toppar_path as global_toppar
from haddock.core.defaults import cns_exec as global_cns_exec
from haddock.core.typing import Any, FilePath, Optional, Union
from haddock.gear.expandable_parameters import populate_mol_parameters_in_module
from haddock.libs.libio import working_directory
from haddock.libs.libutil import sort_numbered_paths
from haddock.modules import BaseHaddockModule


class BaseCNSModule(BaseHaddockModule):
    """
    Operation module for CNS.

    Contains additional functionalities excusive for CNS modules.
    """

    def __init__(self, order: int, path: Path, initial_params: FilePath,
                 cns_script: FilePath) -> None:
        """
        Instantitate CNSModule.

        Parameters
        ----------
        cns_script : str or pathlib.Path
            Path to the main module's cns script.
        """
        super().__init__(order, path, initial_params)

        self.cns_folder_path = Path(Path(cns_script).resolve().parent)
        self.cns_protocol_path = Path(cns_script)
        self.toppar_path = global_toppar
        self.recipe_str = self.cns_protocol_path.read_text()

    def run(self, **params: Any) -> None:
        """Execute the module."""
        log.info(f'Running [{self.name}] module')

        self.update_params(**params)

        # the `mol_*` parameters exist only for CNS jobs.
        if self._num_of_input_molecules:
            populate_mol_parameters_in_module(
                self._params,
                self._num_of_input_molecules,
                self._original_params,
                )

        self.add_parent_to_paths()
        self.envvars = self.default_envvars()

        if self.params['self_contained']:
            self.make_self_contained()

        with working_directory(self.path):
            self._run()

        log.info(f'Module [{self.name}] finished.')

    def default_envvars(self) -> dict[str, str]:
        """Return default env vars updated to `envvars` (if given)."""
        default_envvars = {
            "MODULE": str(self.cns_folder_path),
            "MODDIR": ".",
            "TOPPAR": str(self.toppar_path),
            }

        return default_envvars

    def save_envvars(self, filename: FilePath = "envvars") -> None:
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

    def make_self_contained(self) -> None:
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

        _cns_exec = self.params["cns_exec"] or global_cns_exec
        new_cns = Path(".", Path(_cns_exec).name)
        if not new_cns.exists():
            self.params["cns_exec"] = shutil.copyfile(_cns_exec, new_cns)
            shutil.copystat(_cns_exec, new_cns)
            self.params["cns_exec"] = Path("..", Path(_cns_exec).name)

    def get_ambig_fnames(
            self, prev_ambig_fnames: list[Union[None, FilePath]]
            ) -> Union[list[FilePath], None]:
        """Get the correct ambiguous restraint names.
        
        Parameters
        ----------
        prev_ambig_fnames : list
            list of ambig_fname files encoded in previous models

        Returns
        -------
        ambig_fnames : list or None
            list of ambig_fname files to be used by the CNS module
        """
        ambig_fname: Optional[Path] = self.params["ambig_fname"]
        ambig_fnames = None
        if ambig_fname:
            if ambig_fname.name.endswith("tgz"):
                exp_name = ambig_fname.name.split(".tbl.tgz")[0]
                exp_dir = ambig_fname.parent
                self.log(f"Searching for {exp_name}*tbl files in {exp_dir}")
                path = ambig_fname.parent
                ambig_fnames = list(path.glob(f"{exp_name}*tbl"))
                # abort execution if no files are found
                if len(ambig_fnames) == 0:
                    raise Exception(
                        f"No {exp_name}*tbl files found in {exp_dir}"
                        )
                self.log(f"Found {len(ambig_fnames)} compatible tbl files")
                ambig_fnames = sort_numbered_paths(*ambig_fnames)
        else:
            if self.params["previous_ambig"]:
                # check if there is restraint information in all models
                if None in prev_ambig_fnames:
                    raise Exception("'previous_ambig' option selected but no available restraint information in models")  # noqa: E501
                self.log("Using previously defined restraints")
                ambig_fnames = prev_ambig_fnames.copy()
        return ambig_fnames
