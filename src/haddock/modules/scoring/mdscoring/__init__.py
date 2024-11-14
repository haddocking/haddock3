"""MD scoring module.

This module will perform a short MD simulation on the input models and
score them. No restraints are applied during this step.
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath
from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input, prepare_expected_pdb
from haddock.libs.libontology import PDBFile
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import get_engine
from haddock.modules.scoring import CNSScoringModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(CNSScoringModule):
    """HADDOCK3 module to perform energy minimization scoring."""

    name = RECIPE_PATH.name

    def __init__(
        self, order: int, path: Path, initial_params: FilePath = DEFAULT_CONFIG
    ) -> None:
        cns_script = Path(RECIPE_PATH, "cns", "mdscoring.cns")
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        return

    def _run(self) -> None:
        """Execute module."""
        # Pool of jobs to be executed by the CNS engine
        jobs: list[CNSJob] = []

        try:
            models_to_score = self.previous_io.retrieve_models(individualize=True)
        except Exception as e:
            self.finish_with_error(e)

        self.output_models = []
        for model_num, model in enumerate(models_to_score, start=1):
            scoring_inpyt = prepare_cns_input(
                model_num,
                model,
                self.path,
                self.recipe_str,
                self.params,
                "mdscoring",
                native_segid=True,
                debug=self.params["debug"],
                seed=model.seed if isinstance(model, PDBFile) else None,
            )

            scoring_out = f"mdscoring_{model_num}.out"
            err_fname = f"mdscoring_{model_num}.cnserr"

            # create the expected PDBobject
            expected_pdb = prepare_expected_pdb(model, model_num, ".", "mdscoring")
            # fill the ori_name field of expected_pdb
            expected_pdb.ori_name = model.file_name
            expected_pdb.md5 = model.md5
            expected_pdb.restr_fname = model.restr_fname

            self.output_models.append(expected_pdb)

            job = CNSJob(scoring_inpyt, scoring_out, err_fname, envvars=self.envvars)

            jobs.append(job)

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Get the weights from the defaults
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        # Check for generated output, fail it not all expected files are found
        for pdb in self.output_models:
            if pdb.is_present():
                haddock_model = HaddockModel(pdb.file_name)
                pdb.unw_energies = haddock_model.energies

                haddock_score = haddock_model.calc_haddock_score(**weights)
                pdb.score = haddock_score

        output_fname = "mdscoring.tsv"
        self.log(f"Saving output to {output_fname}")
        self.output(output_fname)
        self.export_io_models(faulty_tolerance=self.params["tolerance"])

        if self.params["per_interface_scoring"]:
            self.per_interface_output(output_fname)
