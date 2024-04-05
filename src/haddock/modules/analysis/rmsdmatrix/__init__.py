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
from pathlib import Path
import os

from haddock import log, RMSD_path
from haddock.core.typing import Any, AtomsDict, FilePath
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libontology import ModuleIO, RMSDFile
from haddock.libs.libparallel import get_index_list
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule
from haddock.modules import get_engine
from haddock.modules.analysis import (
    confirm_resdic_chainid_length,
    get_analysis_exec_mode,
    )
from haddock.modules.analysis.rmsdmatrix.rmsd import (
    RMSDJob,
    rmsd_dispatcher,
    XYZWriter,
    XYZWriterJob
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")
EXEC_PATH = Path(RMSD_path, "src/fast-rmsdmatrix")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with RMSD."""

    name = RECIPE_PATH.name

    def __init__(self,
                 order: int,
                 path: Path,
                 initial_params: FilePath = DEFAULT_CONFIG) -> None:
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if fast-rmsdmatrix is installed and available."""

        if not os.access(EXEC_PATH, mode=os.F_OK):
            raise Exception(f"Required {str(EXEC_PATH)} file does not exist.{os.linesep}"
                            "Old HADDOCK3 installation? Please follow the new installation instructions at https://github.com/haddocking/haddock3/blob/main/docs/INSTALL.md")

        if not os.access(EXEC_PATH, mode=os.X_OK):
            raise Exception(f"Required {str(EXEC_PATH)} file is not executable")

        return

    def _rearrange_output(self, output_name: FilePath, path: FilePath,
                          ncores: int) -> None:
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
        

    def _rearrange_xyz_files(self, output_name: FilePath, path: FilePath,
                             ncores: int) -> None:
        """Combine different xyz outputs in a single file."""
        output_fname = Path(path, output_name)
        self.log(f"rearranging xyz files into {output_fname}")
        # Combine files
        with open(output_fname, 'w') as out_file:
            for core in range(ncores):
                tmp_file = Path(path, "traj_" + str(core) + ".xyz")
                with open(tmp_file) as infile:
                    out_file.write(infile.read())
                log.debug(f"File number {core} written")
                tmp_file.unlink()
        log.info("Completed reconstruction of xyz files.")
        log.info(f"{output_fname} created.")
        

    def check_common_atoms(self, common_keys, coord_keys_lengths):
        # checking the common atoms
        n_atoms = len(common_keys) #common atoms
        max_n_atoms = max(coord_keys_lengths)
        perc = (n_atoms / max_n_atoms) * 100
        if perc == 100.0:
            log.info("All the models share the same atoms.")
        elif perc > self.params["atom_similarity"] and perc < 100.0:
            # if it's between 0.9 and 1, it's likely that the models share the same atoms
            # but still the user may want to see a warning
            log.warning(
                "Not all the atoms are common to all the models."
                f" Common atoms ({n_atoms}) != max_n_atoms {max_n_atoms}. Similarity ({perc:.2f}%) higher than allowed ({self.params['atom_similarity']:.2f}%)."
                )
        else:
            # common keys are less than 90% of the previous keys
            # something is likely wrong
            self.finish_with_error(
                "Input atoms are not the same for all the models."
                f" Common atoms ({n_atoms}) != max_n_atoms {max_n_atoms}. Similarity ({perc:.2f}%) lower than allowed ({self.params['atom_similarity']:.2f}%)."
                " Please check the input ensemble."
                )
        return n_atoms

    def update_params(self, *args: Any, **kwargs: Any) -> None:
        """Update parameters."""
        super().update_params(*args, **kwargs)
        with contextlib.suppress(KeyError):
            self.params.pop("resdic_")

        confirm_resdic_chainid_length(self._params)

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        # Get the models generated in previous step
        models = self.previous_io.retrieve_models(
            individualize=True
            )

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
        # index_list for the jobs with linear scaling
        ncores_linear = parse_ncores(n=self.params['ncores'], njobs=len(models))
        index_list = get_index_list(nmodels, ncores_linear)
        # create prev_keys to check if the keys of the ref_coord_dic
        # are the same for all the models
        common_keys : list[str] = []
        traj_filename = Path("traj.xyz")
        coord_keys_lengths = []

        filter_resdic = {
            key[-1]: value for key, value
            in self.params.items()
            if key.startswith("resdic")
        }
        
        # find common keys
        for mod in models:
            atoms: AtomsDict = get_atoms(mod, self.params["allatoms"])
            
            ref_coord_dic, _ = load_coords(
            mod, atoms, filter_resdic
            )
            coord_keys_lengths.append(len(ref_coord_dic.keys()))
            if common_keys != []:
                common_keys = set(ref_coord_dic.keys()).intersection(common_keys)
            else:
                common_keys = ref_coord_dic.keys()
        
        # check common atoms
        n_atoms = self.check_common_atoms(common_keys, coord_keys_lengths)
        
        xyzwriter_jobs: list[XYZWriterJob] = []
        for core in range(ncores_linear):
            output_name = Path("traj_" + str(core) + ".xyz")
            # init RMSDJobFast
            xyzwriter_obj = XYZWriter(
                model_list=models[index_list[core]:index_list[core + 1]],
                output_name=output_name,
                core=core,
                path=Path("."),
                n_atoms=n_atoms,
                common_keys=common_keys,
                filter_resdic=filter_resdic,
                allatoms=self.params["allatoms"],
                )
            #job_f = output_name
            job = XYZWriterJob(
                xyzwriter_obj,
                )
            xyzwriter_jobs.append(job)
        
        # run jobs
        exec_mode = get_analysis_exec_mode(self.params["mode"])
        Engine = get_engine(exec_mode, self.params)
        engine = Engine(xyzwriter_jobs)
        engine.run()

        self._rearrange_xyz_files(
            traj_filename,
            path=Path("."),
            ncores=ncores_linear
            )
        
        # Calculate the rmsd for each set of models
        rmsd_jobs: list[RMSDJob] = []
        self.log(f"running RmsdFast Jobs with {ncores} cores")
        for core in range(ncores):
            output_name = Path("rmsd_" + str(core) + ".out")
            # init RMSDJobFast
            job = RMSDJob(
                traj_filename,
                output_name,
                EXEC_PATH,
                core,
                npairs[core],
                ref_structs[core],
                mod_structs[core],
                len(models),
                n_atoms,
                )
            rmsd_jobs.append(job)
        
        engine = Engine(rmsd_jobs)
        engine.run()

        rmsd_file_l: list[str] = []
        not_found: list[str] = []
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
        final_output_name = "rmsd.matrix"
        self._rearrange_output(
            final_output_name,
            path=Path("."),
            ncores=ncores
            )
        # Delete the trajectory file
        if traj_filename.exists():
            os.unlink(traj_filename)

        # Sending models to the next step of the workflow
        self.output_models = models
        self.export_io_models()
        # Sending matrix path to the next step of the workflow
        matrix_io = ModuleIO()
        rmsd_matrix_file = RMSDFile(
            final_output_name,
            npairs=tot_npairs
            )
        matrix_io.add(rmsd_matrix_file)
        matrix_io.save(filename="rmsd_matrix.json")
