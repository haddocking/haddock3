"""EM scoring module.

This module performs energy minimization and scoring of the models generated in
the previous step of the workflow. No restraints are applied during this step.
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath
from haddock.libs.libontology import PDBFile
from haddock.modules import get_engine
from haddock.modules.scoring import ScoringModule
from haddock.modules.scoring.prodigy.prodigy import ProdigyJob, CheckInstall


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(ScoringModule):
    """HADDOCK3 module to perform energy minimization scoring."""

    name = RECIPE_PATH.name

    def __init__(
            self,
            order: int,
            path: Path,
            initial_params: FilePath = DEFAULT_CONFIG,
            ) -> None:
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        CheckInstall()
        return

    def _run(self) -> None:
        """Execute module."""
        # Retrieve previous models
        try:
            models_to_score = self.previous_io.retrieve_models(individualize=True)
        except Exception as e:
            self.finish_with_error(e)

        # Pool of Prodigy jobs to be executed
        jobs: list[ProdigyJob] = []
        self.output_models: list[PDBFile] = []

        for mi, model in enumerate(models_to_score):
            self.output_models.append(model)
            prodigy_job = ProdigyJob(
                model,
                self.params,
                index=mi,
                )
            jobs.append(prodigy_job)

        # Run CNS Jobs
        self.log(f"Running {len(jobs)} Prodigy Jobs")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("Prodigy jobs have finished")

        # Check for generated output, fail it not all expected files are found
        for model_score in engine.results:
            # point score stored in engine
            prodigy_score = model_score.score
            model_index = model_score.index
            # Point corresponding model
            self.output_models[model_index].score = prodigy_score

        output_fname = "prodigy.tsv"
        self.log(f"Saving output to {output_fname}")
        self.output(output_fname)

        self.export_io_models()