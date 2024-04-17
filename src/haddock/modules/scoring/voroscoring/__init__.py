"""
"""
from pathlib import Path

from haddock.core.typing import FilePath
from haddock.modules import get_engine
from haddock.modules.scoring import ScoringModule
from haddock.modules.scoring.voroscoring.voroscoring import (
    VoroMQA,
    update_models_with_scores,
    write_models_scores,
    )

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(ScoringModule):
    """."""

    name = RECIPE_PATH.name

    def __init__(self,
                 order: int,
                 path: Path,
                 initial_params: FilePath = DEFAULT_CONFIG) -> None:
        super().__init__(order, path, initial_params=initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        # FIXME ? Check if conda env is accessible
        return

    def _run(self) -> None:
        """Execute module."""
        # Retrieve previous models
        try:
            models_to_score = self.previous_io.retrieve_models(
                individualize=True
                )
        except Exception as e:
            self.finish_with_error(e)

        jobs: list[VoroMQA] = []
        output_fname = f"{RECIPE_PATH.name}_voro.tsv"
        voromqa = VoroMQA(
            models_to_score,
            './',
            self.params,
            output_filepath=output_fname,
            )
        jobs: list[VoroMQA] = [voromqa]

        # Run CNS Jobs
        self.log(f"Running Voro-mqa scoring")
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("Voro-mqa scoring finished!")

        # Update score of output models
        self.output_models, models_scores = update_models_with_scores(
            output_fname,
            models_to_score,
            )
        # Write output file
        scoring_tsv_fpath = f"{RECIPE_PATH.name}.tsv"
        write_models_scores(models_scores, scoring_tsv_fpath)
        # Export to next module
        self.export_io_models(faulty_tolerance=self.params["tolerance"])
