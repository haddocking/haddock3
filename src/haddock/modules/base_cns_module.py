"""Functionalities related to CNS modules."""

import os
import re
import shutil
import tarfile
from pathlib import Path

from haddock import log
from haddock import toppar_path as global_toppar
from haddock.core.defaults import cns_exec as global_cns_exec
from haddock.core.typing import Any, FilePath, Optional, ParamDict, Union
from haddock.gear.expandable_parameters import populate_mol_parameters_in_module
from haddock.libs.libio import working_directory
from haddock.modules import BaseHaddockModule


_CNS_VARIABLE_PATTERN = re.compile(r"\$([A-Za-z0-9_]+)")
_CNS_SPLICED_VARIABLE_PREFIX_PATTERN = re.compile(r"\$([A-Za-z][A-Za-z0-9_]*?)(?=\$)")


class BaseCNSModule(BaseHaddockModule):
    """
    Operation module for CNS.

    Contains additional functionalities excusive for CNS modules.
    """

    # ``ncs_*`` values are consumed by CNS' built-in NCS data structure rather
    # than referenced directly by the module recipe tree.
    CNS_PARAM_INCLUDE_PREFIXES = ("mol_", "fle_", "ncs_")

    def __init__(
        self, order: int, path: Path, initial_params: FilePath, cns_script: FilePath
    ) -> None:
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
        recipe_texts = tuple(
            recipe.read_text(encoding="utf-8")
            for recipe in self.cns_folder_path.rglob("*.cns")
        )
        self._cns_recipe_variables = frozenset(
            variable
            for recipe_text in recipe_texts
            for variable in _CNS_VARIABLE_PATTERN.findall(recipe_text)
        )
        self._cns_spliced_variable_prefixes = tuple(
            sorted(
                {
                    prefix
                    for recipe_text in recipe_texts
                    for prefix in _CNS_SPLICED_VARIABLE_PREFIX_PATTERN.findall(
                        recipe_text
                    )
                },
                key=len,
                reverse=True,
            )
        )

    def run(self, **params: Any) -> None:
        """Execute the module."""
        log.info(f"Running [{self.name}] module")

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

        if self.params["self_contained"]:
            self.make_self_contained()

        with working_directory(self.path):
            self._run()

        log.info(f"Module [{self.name}] finished.")

    def default_envvars(self) -> dict[str, str]:
        """Return default env vars updated to `envvars` (if given)."""
        default_envvars = {
            "MODULE": str(self.cns_folder_path),
            "MODDIR": ".",
            "TOPPAR": str(self.toppar_path),
        }

        return default_envvars

    def cns_params(self, params: Optional[ParamDict] = None) -> ParamDict:
        """Return only parameters that a CNS recipe can consume.

        The inclusion rule is deliberately derived from the module's CNS recipe
        tree rather than maintained as a deny-list of Python orchestration
        settings.  New scheduler or module-control parameters therefore cannot
        silently become part of a CNS job identity.
        """
        source = self.params if params is None else params
        return {
            key: value
            for key, value in source.items()
            if key.rstrip() in self._cns_recipe_variables
            or key.startswith(self.CNS_PARAM_INCLUDE_PREFIXES)
            or self._matches_cns_spliced_variable(key.rstrip())
        }

    def _matches_cns_spliced_variable(self, parameter: str) -> bool:
        """Return whether a parameter can be assembled by a CNS symbol splice."""
        return any(
            re.fullmatch(rf"{re.escape(prefix)}\d+(?:_\d+)*", parameter)
            for prefix in self._cns_spliced_variable_prefixes
        )

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
            shutil.copyfile(_cns_exec, new_cns)
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
        ambig_fnames: Union[list[FilePath], None] = None
        if ambig_fname:
            if ambig_fname.name.endswith("tgz"):
                self.log(f"Searching for tbl files in {ambig_fname}")
                basepath = ambig_fname.parent
                # Open archive to retrieve members
                with tarfile.open(ambig_fname, mode="r:gz") as targz:
                    # Get filenames
                    tbl_files = targz.getmembers()
                # At that stage, they should be already extracted
                # and we simply need to generate the full paths
                ambig_fnames = []
                for tbl in tbl_files:
                    fname = tbl.name
                    # Filter non-".tbl" and hidden files
                    if fname.startswith("."):
                        continue
                    # Build path
                    tbl_fpath = Path(basepath, fname)
                    # Make sure it exists
                    if tbl_fpath.exists() and tbl_fpath.suffix == ".tbl":
                        ambig_fnames.append(tbl_fpath)
                # abort execution if no files are found
                if len(ambig_fnames) == 0:
                    raise Exception(f"No tbl files found in {ambig_fname} !")
                self.log(f"Found {len(ambig_fnames)} tbl files")
                # Sort them to make sure we get the same order every time
                ambig_fnames.sort()
        else:
            if self.params["previous_ambig"]:
                # check if there is restraint information in all models
                if None in prev_ambig_fnames:
                    raise Exception(
                        "'previous_ambig' option selected but no available "
                        "restraint information in models"
                    )
                self.log("Using previously defined restraints")
                ambig_fnames = prev_ambig_fnames.copy()
        return ambig_fnames
