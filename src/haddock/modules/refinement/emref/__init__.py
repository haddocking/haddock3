"""HADDOCK3 module for energy minimization refinement."""
from pathlib import Path

from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input, prepare_expected_pdb
from haddock.libs.libontology import ModuleIO
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import get_engine
from haddock.modules.base_cns_module import BaseCNSModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module energy minimization refinement."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        """."""
        cns_script = Path(RECIPE_PATH, "cns", "emref.cns")
        super().__init__(order, path, initial_params, cns_script=cns_script)

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
        idx = 1
        sampling_factor = self.params["sampling_factor"]
        if sampling_factor > 1:
            self.log(f"sampling_factor={sampling_factor}")
        if sampling_factor == 0:
            self.log("[Warning] sampling_factor cannot be 0, setting it to 1")
            sampling_factor = 1
        if sampling_factor > 100:
            self.log("[Warning] sampling_factor is larger than 100")

        idx = 1
        for model in models_to_refine:
            for _ in range(self.params['sampling_factor']):
                inp_file = prepare_cns_input(
                    idx,
                    model,
                    self.path,
                    self.recipe_str,
                    self.params,
                    "emref",
                    )
                out_file = f"emref_{idx}.out"

                # create the expected PDBobject
                expected_pdb = prepare_expected_pdb(
                    model, idx, ".", "emref"
                    )

                refined_structure_list.append(expected_pdb)

                job = CNSJob(inp_file, out_file, envvars=self.envvars)

                jobs.append(job)

                idx += 1

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Get the weights needed for the CNS module
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        for pdb in refined_structure_list:
            if pdb.is_present():
                haddock_score = HaddockModel(pdb.file_name).calc_haddock_score(
                    **weights
                    )

                pdb.score = haddock_score

        # Save module information
        io = ModuleIO()
        io.add(refined_structure_list, "o")
        faulty = io.check_faulty()
        tolerance = self.params["tolerance"]
        if faulty > tolerance:
            _msg = (
                f"{faulty:.2f}% of output was not generated for this module "
                f"and tolerance was set to {tolerance:.2f}%.")
            self.finish_with_error(_msg)
        io.save()
