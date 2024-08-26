"""
RMSD matrix module.

This module calculates of the RMSD matrix between all the models
generated in the previous step.

As all the pairwise RMSD calculations are independent, the module distributes
them over all the available cores in an optimal way.

Once created, the RMSD matrix is saved in text form in the current `rmsdmatrix`
folder. The path to this file is then shared with the following step of the
workflow by means of the json file `rmsd_matrix.json`.

The module accepts two parameters in input, namely:

* `max_models` (default = 10000)
* `resdic_` : an expandable parameter to specify which residues must be
  considered for the alignment and the RMSD calculation. If there are
  two proteins denoted by chain IDs A and B, then the user can operate
  such selection in the following way inside the configuration file

>>> resdic_A = [1,2,3,4]
>>> resdic_B = [2,3,4]

thus telling the module to consider residues from 1 to 4 of chain A and from 2
to 4 of chain B for the alignment and RMSD calculation.
"""

import contextlib
import os
from pathlib import Path

from haddock import RMSD_path, log
from haddock.core.typing import Any, FilePath
from haddock.libs.libalign import check_common_atoms
from haddock.libs.libio import dump_output_data_to_file
from haddock.libs.libontology import ModuleIO, RMSDFile
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule, get_engine
from haddock.modules.analysis import (
    confirm_resdic_chainid_length,
    get_analysis_exec_mode,
)
from haddock.modules.analysis.rmsdmatrix.rmsd import (
    RMSDJob,
    XYZWriter,
    rmsd_dispatcher,
)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")
EXEC_PATH = Path(RMSD_path, "src/fast-rmsdmatrix")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with RMSD."""

    name = RECIPE_PATH.name

    def __init__(
        self, order: int, path: Path, initial_params: FilePath = DEFAULT_CONFIG
    ) -> None:
        super().__init__(order, path, initial_params)
        self.rmsd_matrix_fname = Path(".", "rmsd.matrix")
        self.traj_filename = Path(".", "traj.xyz")

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if fast-rmsdmatrix is installed and available."""

        if not os.access(EXEC_PATH, mode=os.F_OK):
            raise Exception(
                f"Required {str(EXEC_PATH)} file does not exist.{os.linesep}"
                "Old HADDOCK3 installation? Please follow the new installation instructions at https://github.com/haddocking/haddock3/blob/main/docs/INSTALL.md"
            )

        if not os.access(EXEC_PATH, mode=os.X_OK):
            raise Exception(f"Required {str(EXEC_PATH)} file is not executable")

        return

    def update_params(self, *args: Any, **kwargs: Any) -> None:
        """Update parameters."""
        super().update_params(*args, **kwargs)
        with contextlib.suppress(KeyError):
            self.params.pop("resdic_")

        confirm_resdic_chainid_length(self._params)

    def _run(self) -> None:
        """Execute module."""

        # FIXME: TEMPORARY FIX TO AVOID BATCH MODE
        exec_mode = get_analysis_exec_mode(self.params["mode"])
        if exec_mode == "batch":
            self.finish_with_error("Batch mode not supported for this module")
        Engine = get_engine(exec_mode, self.params)

        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        # Get the models generated in previous step
        models = self.previous_io.retrieve_models(individualize=True)

        nmodels = len(models)
        if nmodels > self.params["max_models"]:
            # too many input models : RMSD matrix would be too big => Abort!
            raise Exception("Too many models for RMSD matrix calculation")

        # index_list for the jobs with linear scaling
        ncores = parse_ncores(n=self.params["ncores"], njobs=len(models))

        filter_resdic = {
            key[-1]: value
            for key, value in self.params.items()
            if key.startswith("resdic")
        }

        # check common atoms
        n_atoms, common_keys = check_common_atoms(
            models,
            filter_resdic,
            self.params["allatoms"],
            self.params["atom_similarity"],
        )

        xyzwriter_jobs: list[XYZWriter] = []

        for model in models:
            xyzwriter_jobs.append(
                XYZWriter(
                    model=model,  # type: ignore
                    n_atoms=n_atoms,
                    common_keys=common_keys,
                    filter_resdic=filter_resdic,
                    allatoms=self.params["allatoms"],
                )
            )

        # run jobs
        engine = Engine(xyzwriter_jobs)
        engine.run()

        # FIXME: This does not work with `HPCSheduler`!
        dump_output_data_to_file(
            input=engine.results, output_fname=self.traj_filename  # type: ignore
        )

        # ============================================================================== #
        # ============================================================================== #
        # TODO: Refactor the code below into a FastRMSDMatrix class
        #  Explanation: `fast-rmsdmatrix` is a standalone executable that
        # calculates a RMSD matrix. It takes as arguments a trajectory file
        # and a "slice" of the matrix to be calculated. The slice is defined by
        # the code below via `rmsd_dispatcher`. Because Python does not do zero-cost
        # abstraction, calculating each line of the matrix in parallel using
        # systems calls to `fast-rmsdmatrix` is NOT the most efficient way to do it.
        # Hence why the logic of the code below should be refactored into a class.
        #
        # This class must work with the Schedulers.
        #
        # ============================================================================== #

        tot_npairs = nmodels * (nmodels - 1) // 2
        log.info(f"total number of pairs {tot_npairs}")
        ncores = parse_ncores(n=self.params["ncores"], njobs=tot_npairs)
        npairs, ref_structs, mod_structs = rmsd_dispatcher(
            nmodels=nmodels, tot_npairs=tot_npairs, ncores=tot_npairs
        )

        # Calculate the rmsd for each set of models
        fast_rmsdmatrix_jobs: list[RMSDJob] = []
        self.log(f"running RmsdFast Jobs with {ncores} cores")
        for core in range(ncores):
            output_name = Path("rmsd_" + str(core) + ".matrix")
            job = RMSDJob(
                self.traj_filename,  # input
                output_name,  # output
                EXEC_PATH,  # executable
                core,
                npairs[core],
                ref_structs[core],
                mod_structs[core],
                len(models),
                n_atoms,
            )
            fast_rmsdmatrix_jobs.append(job)
        # ============================================================================== #
        # ==============================================================================#

        engine = Engine(fast_rmsdmatrix_jobs)
        engine.run()

        # FIXME: This does not work with `HPCSheduler`!
        dump_output_data_to_file(
            input=engine.results, output_fname=self.rmsd_matrix_fname  # type: ignore
        )

        if not self.rmsd_matrix_fname.exists():
            self.finish_with_error("RMSD matrix file was not generated")

        # Sending models to the next step of the workflow
        self.output_models = models  # type: ignore
        self.export_io_models()

        # Sending matrix path to the next step of the workflow
        matrix_io = ModuleIO()
        rmsd_matrix_file = RMSDFile(self.rmsd_matrix_fname, npairs=tot_npairs)
        matrix_io.add(rmsd_matrix_file)
        matrix_io.save(filename="rmsd_matrix.json")
