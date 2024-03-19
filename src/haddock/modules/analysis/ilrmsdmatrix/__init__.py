"""
ilRMSD matrix module.

This module calculates of the interface-ligand RMSD (ilRMSD) matrix between all
 the models generated in the previous step.

As all the pairwise ilRMSD calculations are independent, the module distributes
them over all the available cores in an optimal way.

Once created, the ilRMSD matrix is saved in text form in the current
 `ilrmsdmatrix` folder. The path to this file is then shared with the following
 step of the workflow by means of the json file `rmsd_matrix.json`.

The module accepts one parameters in input, namely:
* `max_models` (default = 10000)

IMPORTANT: the module assumes coherent numbering for all the receptor and
 ligand chains, as no alignment is performed. The user must ensure that the
 numbering is coherent.
"""
import os
from pathlib import Path

import numpy as np

from haddock import log, RMSD_path
from haddock.core.typing import AtomsDict
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libontology import ModuleIO, RMSDFile
from haddock.libs.libparallel import get_index_list
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule, read_from_yaml_config
from haddock.modules import get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.ilrmsdmatrix.ilrmsd import (
    Contact,
    ContactJob,
    )
from haddock.modules.analysis.rmsdmatrix import rmsd_dispatcher, RMSDJob



RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")
EXEC_PATH = Path(RMSD_path, "src/fast-rmsdmatrix")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with RMSD."""

    name = RECIPE_PATH.name
    
    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if fast-rmsdmatrix is installed and available."""

        if not os.access(EXEC_PATH, mode=os.F_OK):
            raise Exception(f"Required {str(EXEC_PATH)} file does not exist.{os.linesep}"
                            "Old HADDOCK3 installation? Please follow the new installation instructions at https://github.com/haddocking/haddock3/blob/main/docs/INSTALL.md")

        if not os.access(EXEC_PATH, mode=os.X_OK):
            raise Exception(f"Required {str(EXEC_PATH)} file is not executable")

        return

    @staticmethod
    def _rearrange_output(output_name, path, ncores):
        """Combine different ilrmsd outputs in a single file."""
        output_fname = Path(path, output_name)
        log.info(f"rearranging output files into {output_fname}")
        # Combine files
        with open(output_fname, 'w') as out_file:
            for core in range(ncores):
                tmp_file = Path(path, "ilrmsd_" + str(core) + ".matrix")
                with open(tmp_file) as infile:
                    out_file.write(infile.read())
                log.debug(f"File number {core} written")
                tmp_file.unlink()
        log.info("Completed reconstruction of rmsd files.")
        log.info(f"{output_fname} created.")

    @staticmethod
    def _rearrange_contact_output(output_name, path, ncores):
        """Combine different contact outputs in a single file.
        
        Parameters
        ----------
        output_name : str
            Name of the output file.
        path : str
            Path to the output file.
        ncores : int
            Number of cores used in the calculation.
        
        Returns
        -------
        res_resdic_npu : dict
            Dictionary of the unique residues in the interfaces.
        """
        output_fname = Path(path, output_name)
        log.info(f"rearranging contact files into {output_fname}")
        # Combine residues
        res_resdic = {}
        for core in range(ncores):
            tmp_file = Path(path, "interface_contacts_" + str(core) + ".con")
            for ln in tmp_file.read_text().split(os.linesep)[:-1]:
                ch = ln.split()[0]
                resids = [int(el) for el in ln.split()[1:]]
                if ch not in res_resdic:
                    res_resdic[ch] = [resids]
                else:
                    res_resdic[ch].append(resids)
            tmp_file.unlink()
        
        # Unique residues
        npu_resdic = {}
        for ch in res_resdic.keys():
            resids = np.concatenate(res_resdic[ch])
            resids_npu = np.unique(resids)
            npu_resdic[ch] = resids_npu
        log.info(f"Overall interface residues: {npu_resdic}")

        # Write to file
        with open(output_fname, 'w') as out_file:
            for chain in npu_resdic.keys():
                out_file.write(f"{chain} ")
                out_file.write(" ".join([str(el) for el in npu_resdic[chain]]))
                out_file.write(os.linesep)

        log.info("Completed reconstruction of contact files.")
        log.info(f"{output_fname} created.")
        return npu_resdic

    def _run(self):
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
        ncores = parse_ncores(n=self.params['ncores'], njobs=nmodels)
        
        if nmodels > self.params['max_models']:
            # too many input models : ilRMSD matrix would be too big => Abort!
            raise Exception(f"Too many models ({nmodels} > {self.params['max_models']}) for ilRMSD matrix calculation")

        # Find the common residues making contacts for the receptor and ligand.
        index_list = get_index_list(nmodels, ncores)
        # contact jobs
        contact_jobs = []
        for core in range(ncores):
            output_name = "interface_contacts_" + str(core) + ".con"
            contact_obj = Contact(
                model_list=models[index_list[core]:index_list[core + 1]],
                output_name=output_name,
                core=core,
                path=Path("."),
                params=self.params,
                )
            # running now the ContactJob
            job_f = Path(output_name)
            # init RMSDJob
            job = ContactJob(
                job_f,
                self.params,
                contact_obj,
                )
            contact_jobs.append(job)

        exec_mode = get_analysis_exec_mode(self.params["mode"])

        Engine = get_engine(exec_mode, self.params)
        contact_engine = Engine(contact_jobs)
        contact_engine.run()

        # check if all the jobs have been completed
        contact_file_l = []
        not_found = []
        for job in contact_jobs:
            if not job.output.exists():
                not_found.append(job.output.name)
                wrn = f'Contacts were not calculated for {job.output.name}'
                log.warning(wrn)
            else:
                contact_file_l.append(str(job.output))
        
        if not_found:
            # Not all contacts were calculated, cannot proceed
            self.finish_with_error("Several contact files were not generated:"
                                   f" {not_found}")
        
        # Post-processing : single file
        output_name = "receptor_contacts.con"
        res_resdic = self._rearrange_contact_output(
            output_name,
            path=contact_obj.path,
            ncores=ncores
            )
        
        tot_npairs = nmodels * (nmodels - 1) // 2
        ncores = parse_ncores(n=self.params['ncores'], njobs=tot_npairs)
        log.info(f"total number of pairs {tot_npairs}")
        npairs, ref_structs, mod_structs = rmsd_dispatcher(
            nmodels,
            tot_npairs,
            ncores)
        
        rec_traj_filename = Path("traj_rec.xyz")
        lig_traj_filename = Path("traj_lig.xyz")
        prev_keys = []
        with open(rec_traj_filename, "w") as rec_traj_xyz:
            with open(lig_traj_filename, "w") as lig_traj_xyz:
                for mod in models:
                    atoms: AtomsDict = get_atoms(mod)
                    ref_coord_dic, _ = load_coords(
                    mod, atoms, res_resdic
                    )
                    if prev_keys != []:
                        if ref_coord_dic.keys() != prev_keys:
                            self.finish_with_error(
                                "The keys of the ref_coord_dic are not the "
                                "same for all the models. Please check the "
                                "input models."
                                )
                    # write receptor coords
                    rec_coords = {k:ref_coord_dic[k] for k in ref_coord_dic if k[0] == self.params["receptor_chain"]}
                    lig_coords = {k:ref_coord_dic[k] for k in ref_coord_dic if k[0] in self.params["ligand_chains"]}
                    n_rec_atoms = len(rec_coords)
                    n_lig_atoms = len(lig_coords)
                    rec_traj_xyz.write(f"{n_rec_atoms}{os.linesep}{os.linesep}")
                    lig_traj_xyz.write(f"{n_lig_atoms}{os.linesep}{os.linesep}")
                    for v in rec_coords.values():
                        rec_traj_xyz.write(f"x {v[0]} {v[1]} {v[2]}{os.linesep}")
                    for v in lig_coords.values():
                        lig_traj_xyz.write(f"x {v[0]} {v[1]} {v[2]}{os.linesep}")
                    
                    prev_keys = ref_coord_dic.keys()
        
        # Calculate the rmsd for each set of models
        ilrmsd_jobs: list[RMSDJob] = []
        self.log(f"running RmsdFast Jobs with {ncores} cores")
        for core in range(ncores):
            output_name = Path("ilrmsd_" + str(core) + ".out")
            job_f = Path(output_name)
            # init RMSDJob
            job = RMSDJob(
                rec_traj_filename,
                output_name,
                EXEC_PATH,
                core,
                npairs[core],
                ref_structs[core],
                mod_structs[core],
                len(models),
                n_rec_atoms,
                lig_traj_filename,
                n_lig_atoms,
                )
            ilrmsd_jobs.append(job)

        Engine = get_engine(exec_mode, self.params)
        ilrmsd_engine = Engine(ilrmsd_jobs)
        ilrmsd_engine.run()

        ilrmsd_file_l = []
        not_found = []
        for j in ilrmsd_jobs:
            if not j.output.exists():
                # NOTE: If there is no output, most likely the RMSD calculation
                # timed out
                not_found.append(j.output.name)
                wrn = f'ilRMSD results were not calculated for {j.output.name}'
                log.warning(wrn)
            else:
                ilrmsd_file_l.append(str(j.output))

        if not_found:
            # Not all distances were calculated, cannot create the full matrix
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        # Post-processing : single file
        output_name = "ilrmsd.matrix"
        self._rearrange_output(
            output_name,
            path=Path("."),
            ncores=ncores
            )
        # Delete the trajectory files
        if rec_traj_filename.exists():
            os.unlink(rec_traj_filename)
        if lig_traj_filename.exists():
            os.unlink(lig_traj_filename)

        # Sending models to the next step of the workflow
        self.output_models = models
        self.export_io_models()
        # Sending matrix path to the next step of the workflow
        matrix_io = ModuleIO()
        ilrmsd_matrix_file = RMSDFile(
            output_name,
            npairs=tot_npairs
            )
        matrix_io.add(ilrmsd_matrix_file)
        matrix_io.save(filename="rmsd_matrix.json")
