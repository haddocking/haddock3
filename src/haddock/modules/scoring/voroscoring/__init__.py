"""Voro scoring module.

This module performs scoring of input pdb models using ftdmp voro-mqa-all tool.
For more information, please check: https://github.com/kliment-olechnovic/ftdmp

It is a third party module, and requires the appropriate set up and intallation
for it to run without issue.
"""
from os import linesep
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.modules import get_engine
from haddock.modules.scoring import ScoringModule
from haddock.modules.scoring.voroscoring.voroscoring import (
    VoroMQA,
    update_models_with_scores,
    )

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(ScoringModule):
    """."""

    name = RECIPE_PATH.name

    def __init__(
            self,
            order: int,
            path: Path,
            *ignore: Any,
            init_params: FilePath = DEFAULT_CONFIG,
            **everything: Any,
            ) -> None:
        """Initialize class."""
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        # FIXME ? Check :
        # - if conda env is accessible
        # - if ftdmp is accessible
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

        # Initiate VoroMQA object
        output_fname = Path("voro_mqa_all.tsv")
        voromqa = VoroMQA(
            models_to_score,
            './',
            self.params,
            output=output_fname,
            )

        # Launch machinery
        jobs: list[VoroMQA] = [voromqa]
        # Run Job(s)
        self.log("Running Voro-mqa scoring")
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("Voro-mqa scoring finished!")

        # Update score of output models
        try:
            self.output_models = update_models_with_scores(
                output_fname,
                models_to_score,
                metric=self.params["metric"],
                )
        except ValueError as e:
            self.finish_with_error(e)

        # Write output file
        scoring_tsv_fpath = f"{RECIPE_PATH.name}.tsv"
        self.output(
            scoring_tsv_fpath,
            header_comments=f"# Note that negative of the value are reported in the case of non-energetical predictions{linesep}",  # noqa : E501
            )
        # Export to next module
        self.export_io_models()