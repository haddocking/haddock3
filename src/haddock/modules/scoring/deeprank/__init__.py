"""Deeprank scoring module.

This module scores complexes using deeprank-gnn-esm, a graph neural network
that predicts the fraction of native contacts (fnat) between two chains of
an interface.

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/scoring.html>`_
"""

from pathlib import Path
from haddock.core.typing import FilePath, Any
from haddock.core.defaults import MODULE_DEFAULT_YAML

from haddock.libs.libontology import PDBFile
from haddock.libs.libutil import parse_ncores
from haddock.modules.scoring import ScoringModule

from haddock.modules.scoring.deeprank.deeprank import (
    DeeprankWrapper,
    deeprank_is_available,
)

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(ScoringModule):
    """HADDOCK3 module to perform deeprank-gnn-esm scoring."""

    name = RECIPE_PATH.name

    def __init__(
        self,
        order: int,
        path: Path,
        *ignore: Any,
        init_params: FilePath = DEFAULT_CONFIG,
        **everything: Any,
    ) -> None:
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if the module is ready to use."""

        if not deeprank_is_available():
            raise Exception(
                "You are trying to use the `deeprank` module but it is not available, please check the installation instructions"
            )

        return

    def _run(self) -> None:
        """Execute module."""
        # individualize=True: deeprank scores each model separately
        #  retrieve_models must return flat PDBFile
        models_to_use = self.previous_io.retrieve_models(individualize=True)
        model_paths = [Path(m.path, m.file_name) for m in models_to_use]

        # NOTE: deeprank has its own logic of parallelization mechanism
        #  so here we DO NOT use haddock's engine and we let deeprank the execution.
        #  Because of that we need `parse_ncores` explicitly
        ncores = parse_ncores(self.params["ncores"])

        deeprank_wrapper = DeeprankWrapper(
            models=model_paths,
            ncores=ncores,
            chain_i=self.params["chain_i"],
            chain_j=self.params["chain_j"],
            path=".",
        )

        deeprank_wrapper.run()
        result_dic: dict[str, float] = deeprank_wrapper.retrieve_scores()

        # Add the score obtained by deeprank back to the models
        for model, model_path in zip(models_to_use, model_paths):
            model.score = result_dic[str(model_path)]

        # Pass the models ahead
        self.output_models: list[PDBFile] = models_to_use

        # Generate a tsv file containing the computed scores
        output_fname = "deeprank.tsv"
        self.log(f"Saving output to {output_fname}")
        self.output(output_fname)

        self.export_io_models()
