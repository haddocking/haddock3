"""HADDOCK3 rigid-body docking module."""
from pathlib import Path

from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input
from haddock.libs.libhpc import HPCScheduler
from haddock.libs.libontology import ModuleIO, PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for rigid body sampling."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = Path(RECIPE_PATH, "cns", "rigidbody.cns")
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
            if self.params["crossdock"]:
                self.log("crossdock=true")
            models_to_dock = self.previous_io.retrieve_models(
                crossdock=self.params["crossdock"]
                )
        except Exception as e:
            self.finish_with_error(e)

        # How many times each combination should be sampled,
        #  cannot be smaller than 1
        sampling_factor = int(self.params["sampling"] / len(models_to_dock))
        if sampling_factor < 1:
            self.finish_with_error(
                "Sampling is smaller than the number"
                " of model combinations "
                f"#model_combinations={len(models_to_dock)},"
                f' sampling={self.params["sampling"]}.'
                )

        # Prepare the jobs
        idx = 1
        structure_list = []
        self.log("Preparing jobs...")
        for combination in models_to_dock:

            for _i in range(sampling_factor):
                inp_file = prepare_cns_input(
                    idx,
                    combination,
                    self.path,
                    self.recipe_str,
                    self.params,
                    "rigidbody",
                    ambig_fname=self.params["ambig_fname"],
                    )

                log_fname = Path(self.path, f"rigidbody_{idx}.out")
                output_pdb_fname = Path(self.path, f"rigidbody_{idx}.pdb")

                # Create a model for the expected output
                model = PDBFile(output_pdb_fname, path=self.path)
                model.topology = [e.topology for e in combination]
                structure_list.append(model)

                job = CNSJob(
                    inp_file,
                    log_fname,
                    cns_folder=self.cns_folder_path,
                    modpath=self.path,
                    config_path=self.params["config_path"],
                    cns_exec=self.params["cns_exec"],
                    )
                jobs.append(job)

                idx += 1

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Get the weights according to CNS parameters
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        not_present = []
        for model in structure_list:
            if not model.is_present():
                not_present.append(model.name)
            else:
                # Score the model
                haddock_score = HaddockModel(
                    model.full_name).calc_haddock_score(
                    **weights
                    )

                model.score = haddock_score

        # Check for generated output
        if len(not_present) == len(structure_list):
            # fail if not all expected files are found
            self.finish_with_error("No models were generated.")

        if not_present:
            # fail if any expected files are found
            self.finish_with_error(
                f"Several models were not generated" f" {not_present}"
                )

        # Save module information
        io = ModuleIO()
        io.add(structure_list, "o")
        io.save(self.path)
