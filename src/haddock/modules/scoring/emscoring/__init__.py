"""EM scoring module.

This module performs energy minimization and scoring of the models generated in
the previous step of the workflow. No restraints are applied during this step.
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath
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
        cns_script = Path(RECIPE_PATH, "cns", f"{self.name}.cns")
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

        # Here we pop the parameter as not supported by CNS and only used
        # at the python level for downstream analysis
        interface_combinations = self.extract_interface_combinations()

        # Itereate over models to prepare CNS inputs
        self.output_models = []
        for model_num, model in enumerate(models_to_score, start=1):
            scoring_input = prepare_cns_input(
                model_num,
                model,
                self.path,
                self.recipe_str,
                self.params,
                self.name,
                native_segid=True,
                debug=self.params["debug"],
                seed=model.seed if isinstance(model, PDBFile) else None,
            )

            scoring_out = f"{self.name}_{model_num}.out"
            err_fname = f"{self.name}_{model_num}.cnserr"

            # create the expected PDBobject
            expected_pdb = prepare_expected_pdb(model, model_num, ".", self.name)
            # fill the ori_name field of expected_pdb
            expected_pdb.ori_name = model.file_name
            expected_pdb.md5 = model.md5
            expected_pdb.restr_fname = model.restr_fname

            self.output_models.append(expected_pdb)

            job = CNSJob(scoring_input, scoring_out, err_fname, envvars=self.envvars)

            jobs.append(job)

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Update the score attributes for each output pdb
        output_haddock_models = self.update_pdb_scores(interface_combinations)
    
        # Set output filename
        output_fname = f"{self.name}.tsv"
        # Process per-interface analysis
        self.per_interface_output(output_fname, output_haddock_models)
        # Generate output
        self.log(f"Saving output to {output_fname}")
        self.output(output_fname)
        # Export models to next module
        self.export_io_models(faulty_tolerance=self.params["tolerance"])
