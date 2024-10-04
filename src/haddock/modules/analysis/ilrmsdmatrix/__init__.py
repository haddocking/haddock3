"""
ilRMSD matrix module.

This module calculates of the interface-ligand RMSD (ilRMSD) matrix between all
the models generated in the previous step.

As all the pairwise ilRMSD calculations are independent, the module distributes
them over all the available cores in an optimal way.

Once created, the ilRMSD matrix is saved in text form in the current
`ilrmsdmatrix` folder. The path to this file is then shared with the following
step of the workflow by means of the json file `rmsd_matrix.json`.

IMPORTANT: the module assumes coherent numbering for all the receptor and ligand
chains, as no alignment is performed. The user must ensure that the numbering
is coherent.
"""

import os
from pathlib import Path

import numpy as np

from haddock import RMSD_path, log
from haddock.core.defaults import FAST_RMSDMATRIX_EXEC, MODULE_DEFAULT_YAML
from haddock.libs.libalign import (
    check_chains,
    check_common_atoms,
    get_atoms,
    load_coords,
    rearrange_xyz_files,
    )
from haddock.libs.libontology import ModuleIO, RMSDFile
from haddock.libs.libparallel import get_index_list
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule, get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.ilrmsdmatrix.ilrmsd import Contact, ContactJob
from haddock.modules.analysis.rmsdmatrix import RMSDJob, rmsd_dispatcher
from haddock.modules.analysis.rmsdmatrix.rmsd import XYZWriter, XYZWriterJob


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)
EXEC_PATH = FAST_RMSDMATRIX_EXEC


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with RMSD."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG, **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if fast-rmsdmatrix is installed and available."""
        if not os.access(EXEC_PATH, mode=os.F_OK):
            raise Exception(
                f"Required {str(EXEC_PATH)} file does not exist.{os.linesep}"
                "Old HADDOCK3 installation? Please follow the new installation instructions at https://github.com/haddocking/haddock3/blob/main/docs/INSTALL.md"  # noqa : E501
            )

        if not os.access(EXEC_PATH, mode=os.X_OK):
            raise Exception(f"Required {str(EXEC_PATH)} file is not executable")

        return

    @staticmethod
    def _rearrange_output(output_name, path, ncores):
        """Combine different ilrmsd outputs in a single file."""
        output_fname = Path(path, output_name)
        log.info(f"rearranging output files into {output_fname}")
        # Combine files
        with open(output_fname, "w") as out_file:
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
        with open(output_fname, "w") as out_file:
            for chain in npu_resdic.keys():
                out_file.write(f"{chain} ")
                out_file.write(" ".join([str(el) for el in npu_resdic[chain]]))
                out_file.write(os.linesep)

        log.info("Completed reconstruction of contact files.")
        log.info(f"{output_fname} created.")
        return npu_resdic

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        # Get the models generated in previous step
        models = self.previous_io.retrieve_models(individualize=True)

        # Parallelisation : optimal dispatching of models
        nmodels = len(models)
        ncores = parse_ncores(n=self.params["ncores"], njobs=nmodels)

        if nmodels > self.params["max_models"]:
            # too many input models : ilRMSD matrix would be too big => Abort!
            raise Exception(
                f"Too many models ({nmodels} > {self.params['max_models']}) "
                "for ilRMSD matrix calculation"
            )

        # find the existing chains
        model_atoms = get_atoms(models[0])
        mod_coord_dic, _ = load_coords(
            models[0],
            model_atoms,
        )
        obs_chains = np.unique([el[0] for el in mod_coord_dic.keys()])
        obs_chains = list(obs_chains)
        log.info(f"Observed chains: {obs_chains}")
        # assigning the chains to the receptor and ligand
        r_chain, l_chains = check_chains(
            obs_chains, self.params["receptor_chain"], self.params["ligand_chains"]
        )
        log.info(f"Receptor chain: {r_chain}")
        log.info(f"Ligand chains: {l_chains}")
        self.params["receptor_chain"] = r_chain
        self.params["ligand_chains"] = l_chains

        # Find the common residues making contacts for the receptor and ligand.
        index_list = get_index_list(nmodels, ncores)
        # contact jobs
        contact_jobs = []
        for core in range(ncores):
            output_name = "interface_contacts_" + str(core) + ".con"
            contact_obj = Contact(
                model_list=models[index_list[core] : index_list[core + 1]],
                output_name=output_name,
                core=core,
                path=Path("."),
                contact_distance_cutoff=self.params["contact_distance_cutoff"],
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
                wrn = f"Contacts were not calculated for {job.output.name}"
                log.warning(wrn)
            else:
                contact_file_l.append(str(job.output))

        if not_found:
            # Not all contacts were calculated, cannot proceed
            self.finish_with_error(
                "Several contact files were not generated:" f" {not_found}"
            )

        # Post-processing : single file
        output_name = "receptor_contacts.con"
        res_resdic = self._rearrange_contact_output(
            output_name, path=contact_obj.path, ncores=ncores
        )

        # if the receptor chain in res_resdic is empty, then the receptor has made no contacts and
        # the ilrmsd matrix cannot be calculated. This probably means that single chains structures
        # reached this step or that something went (very) wrong in the docking process.
        if res_resdic[r_chain].size == 0:
            _msg = f"No contacts found for receptor chain {r_chain}. Impossible to calculate ilRMSD matrix."
            _msg += " Please check your input and make sure that there are at least two chains in contact."
            self.finish_with_error(_msg)

        rec_traj_filename = Path("traj_rec.xyz")
        lig_traj_filename = Path("traj_lig.xyz")

        res_resdic_rec = {k: res_resdic[k] for k in res_resdic if k[0] == r_chain}
        # ligand_chains is a list of chains
        res_resdic_lig = {k: res_resdic[k] for k in l_chains}

        log.info(
            f"Check common atoms for receptor (chain {list(res_resdic_rec.keys())})"
        )
        n_atoms_rec, common_keys_rec = check_common_atoms(
            models,
            res_resdic_rec,
            self.params["allatoms"],
            self.params["atom_similarity"],
        )

        log.info(
            f"Check common atoms for ligand (chains {list(res_resdic_lig.keys())})"
        )
        n_atoms_lig, common_keys_lig = check_common_atoms(
            models,
            res_resdic_lig,
            self.params["allatoms"],
            self.params["atom_similarity"],
        )

        xyzwriter_jobs: list[XYZWriterJob] = []
        for core in range(ncores):
            output_name_rec = Path("traj_rec_" + str(core) + ".xyz")
            # init XYZWriter
            xyzwriter_obj_rec = XYZWriter(
                model_list=models[index_list[core] : index_list[core + 1]],
                output_name=output_name_rec,
                core=core,
                n_atoms=n_atoms_rec,
                common_keys=common_keys_rec,
                filter_resdic=res_resdic_rec,
                allatoms=self.params["allatoms"],
            )
            # job_rec
            job_rec = XYZWriterJob(
                xyzwriter_obj_rec,
            )

            xyzwriter_jobs.append(job_rec)

            output_name_lig = Path("traj_lig_" + str(core) + ".xyz")
            # init XYZWriter
            xyzwriter_obj_lig = XYZWriter(
                model_list=models[index_list[core] : index_list[core + 1]],
                output_name=output_name_lig,
                core=core,
                n_atoms=n_atoms_lig,
                common_keys=common_keys_lig,
                filter_resdic=res_resdic_lig,
                allatoms=self.params["allatoms"],
            )
            # job_lig
            job_lig = XYZWriterJob(
                xyzwriter_obj_lig,
            )
            xyzwriter_jobs.append(job_lig)

        # run jobs
        engine = Engine(xyzwriter_jobs)
        engine.run()

        rearrange_xyz_files(rec_traj_filename, path=Path("."), ncores=ncores)
        rearrange_xyz_files(lig_traj_filename, path=Path("."), ncores=ncores)

        # Parallelisation : optimal dispatching of models
        tot_npairs = nmodels * (nmodels - 1) // 2
        ncores = parse_ncores(n=self.params["ncores"], njobs=tot_npairs)
        log.info(f"total number of pairs {tot_npairs}")
        npairs, ref_structs, mod_structs = rmsd_dispatcher(nmodels, tot_npairs, ncores)

        # Calculate the rmsd for each set of models
        ilrmsd_jobs: list[RMSDJob] = []
        self.log(f"running RmsdFast Jobs with {ncores} cores")
        for core in range(ncores):
            output_name = f"ilrmsd_{core:d}.out"
            job_f = Path(output_name)
            # init RMSDJob
            job = RMSDJob(
                rec_traj_filename,
                job_f,
                EXEC_PATH,
                core,
                npairs[core],
                ref_structs[core],
                mod_structs[core],
                len(models),
                n_atoms_rec,
                lig_traj_filename,
                n_atoms_lig,
            )
            ilrmsd_jobs.append(job)

        ilrmsd_engine = Engine(ilrmsd_jobs)
        ilrmsd_engine.run()

        ilrmsd_file_l = []
        not_found = []
        for j in ilrmsd_jobs:
            if not j.output.exists():
                # NOTE: If there is no output, most likely the RMSD calculation
                # timed out
                not_found.append(j.output.name)
                wrn = f"ilRMSD results were not calculated for {j.output.name}"
                log.warning(wrn)
            else:
                ilrmsd_file_l.append(str(j.output))

        if not_found:
            # Not all distances were calculated, cannot create the full matrix
            self.finish_with_error("Several files were not generated:" f" {not_found}")

        # Post-processing : single file
        output_name = "ilrmsd.matrix"
        self._rearrange_output(output_name, path=Path("."), ncores=ncores)
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
        ilrmsd_matrix_file = RMSDFile(output_name, npairs=tot_npairs)
        matrix_io.add(ilrmsd_matrix_file)
        matrix_io.save(filename="rmsd_matrix.json")
