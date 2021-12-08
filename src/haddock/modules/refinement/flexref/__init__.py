"""HADDOCK3 module for flexible refinement."""
from pathlib import Path

from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input, prepare_expected_pdb
from haddock.libs.libontology import ModuleIO
from haddock.libs.libparallel import Scheduler
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for flexible refinement."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = Path(RECIPE_PATH, "cns", "flexref.cns")
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        """Confirm module is installed."""
        return

    def _run(self):
        """Execute module."""
        # Pool of jobs to be executed by the CNS engine
        jobs = []

        # Get the models generated in previous step
        try:
            models_to_refine = self.previous_io.retrieve_models()
        except Exception as e:
            self.finish_with_error(e)

        refined_structure_list = []
        for idx, model in enumerate(models_to_refine, start=1):
            inp_file = prepare_cns_input(
                idx,
                model,
                self.path,
                self.recipe_str,
                self.params,
                "flexref",
                ambig_fname=self.params["ambig_fname"],
                )

            out_file = Path(self.path, f"flexref_{idx}.out")

            # create the expected PDBobject
            expected_pdb = prepare_expected_pdb(
                model, idx, self.path, "flexref"
                )
            refined_structure_list.append(expected_pdb)

            job = CNSJob(
                inp_file,
                out_file,
                cns_folder=self.cns_folder_path,
                modpath=self.path,
                config_path=self.params["config_path"],
                cns_exec=self.params["cns_exec"],
                )

            jobs.append(job)

        # Run CNS engine
        self.log(f"Running CNS engine with {len(jobs)} jobs")
        engine = Scheduler(jobs, ncores=self.params["ncores"])
        engine.run()
        self.log("CNS engine has finished")

        # Get the weights from the defaults
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        expected = []
        not_found = []
        for pdb in refined_structure_list:
            if not pdb.is_present():
                not_found.append(pdb.file_name)
            else:
                haddock_score = HaddockModel(pdb.full_name).calc_haddock_score(
                    **weights
                    )

                pdb.score = haddock_score
                expected.append(pdb)

        if not_found:
            # fail if any expected files are found
            self.finish_with_error(
                "Several files were not generated:" f" {not_found}"
                )

        # Save module information
        io = ModuleIO()
        io.add(expected, "o")
        io.save(self.path)
