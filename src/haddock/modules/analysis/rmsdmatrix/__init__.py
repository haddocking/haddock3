"""HADDOCK3 RMSD matrix module."""
from pathlib import Path

from haddock import log
from haddock.libs.libontology import ModuleIO, RMSDFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.rmsdmatrix.rmsd import (
    RMSD,
    RMSDJob,
    rmsd_dispatcher,
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with RMSD."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if contact executable is compiled."""
        return

    def _rearrange_output(self, output_name, path, ncores):
        """Combine different rmsd outputs in a single file."""
        output_fname = Path(path, output_name)
        self.log(f"rearranging output files into {output_fname}")
        # Combine files
        with open(output_fname, 'w') as out_file:
            for core in range(ncores):
                tmp_file = Path(path, "rmsd_" + str(core) + ".matrix")
                with open(tmp_file) as infile:
                    out_file.write(infile.read())
                log.debug(f"File number {core} written")
                tmp_file.unlink()
        log.info("Completed reconstruction of rmsd files.")
        log.info(f"{output_fname} created.")

    def _run(self):
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        # Get the models generated in previous step
        models = self.previous_io.retrieve_models()
        
        # Parallelisation : optimal dispatching of models
        nmodels = len(models)
        if nmodels > self.params["max_models"]:
            # too many input models : RMSD matrix would be too big => Abort!
            raise Exception("Too many models for RMSD matrix calculation")
        tot_npairs = nmodels * (nmodels - 1) // 2
        log.info(f"total number of pairs {tot_npairs}")
        ncores = parse_ncores(n=self.params['ncores'], njobs=tot_npairs)
        npairs, ref_structs, mod_structs = rmsd_dispatcher(
            nmodels,
            tot_npairs,
            ncores)

        # Calculate the rmsd for each set of models
        rmsd_jobs = []
        self.log(f"running Rmsd Jobs with {ncores} cores")
        for core in range(ncores):
            output_name = "rmsd_" + str(core) + ".matrix"
            rmsd_obj = RMSD(
                models,
                core,
                npairs[core],
                ref_structs[core],
                mod_structs[core],
                output_name,
                path=Path("."),
                params=self.params
                )
            job_f = Path(output_name)
            # init RMSDJob
            job = RMSDJob(
                job_f,
                self.params,
                rmsd_obj
                )
            rmsd_jobs.append(job)

        rmsd_engine = Scheduler(rmsd_jobs, ncores=ncores)
        rmsd_engine.run()

        rmsd_file_l = []
        not_found = []
        for job in rmsd_jobs:
            if not job.output.exists():
                # NOTE: If there is no output, most likely the RMSD calculation
                # timed out
                not_found.append(job.output.name)
                wrn = f'Rmsd results were not calculated for {job.output.name}'
                log.warning(wrn)
            else:
                rmsd_file_l.append(str(job.output))

        if not_found:
            # Not all distances were calculated, cannot create the full matrix
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        # Post-processing : single file
        output_name = "rmsd.matrix"
        self._rearrange_output(
            output_name,
            path=rmsd_obj.path,
            ncores=ncores
            )

        # Sending models to the next step of the workflow
        self.output_models = models
        self.export_output_models()
        # Sending matrix path to the next step of the workflow
        matrix_io = ModuleIO()
        rmsd_matrix_file = RMSDFile(
            output_name,
            npairs=tot_npairs
            )
        matrix_io.add(rmsd_matrix_file)
        matrix_io.save(filename="rmsd_matrix.json")
