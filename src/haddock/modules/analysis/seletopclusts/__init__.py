"""Select models from the top clusters.

This module selects a number of models from a number of clusters. The
selection is based on the score of the models within the clusters.

In the standard HADDOCK analysis, the top 4 models of the top 10 clusters
are shown. In case seletopclusts is run after a sampling module, we can
keep a few models from all the clusters to have more diversity at the
refinement stage(s).
"""

from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.seletopclusts.seletopclusts import (
    select_top_clusts_models,
    write_selected_models,
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """Haddock Module for 'seletopclusts'."""

    name = RECIPE_PATH.name

    def __init__(self,
                 order: int,
                 path: Path,
                 *ignore: Any,
                 init_params: FilePath = DEFAULT_CONFIG,
                 **everything: Any) -> None:
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if module is installed."""
        return

    def _run(self) -> None:
        """Execute the module's protocol."""
        # Check parameters validity
        if self.params["top_models"] <= 0:
            _msg = "top_models must be either > 0 or nan."
            self.finish_with_error(_msg)

        if not isinstance(self.params["top_clusters"], int):
            _msg = "top_clusters must be an integer."
            self.finish_with_error(_msg)

        # Retrieve list of previous models
        models_to_select = self.previous_io.retrieve_models()

        # Check if cluster info is accessible
        if any([mdl.clt_rank is None for mdl in models_to_select]):
            _msg = (
                "Impossible to obtain cluster information. Please consider "
                "running a clustering method prior to this module."
                )
            self.finish_with_error(_msg)

        # Make model selection
        selected_models, _notes = select_top_clusts_models(
            self.params["sortby"],
            models_to_select,
            self.params["top_clusters"],
            self.params["top_models"],
            )
        # Log notes
        for note in _notes:
            log.info(note)

        # dump the models to disk and change their attributes
        renamed_models = write_selected_models(
            "seletopclusts.txt",
            selected_models,
            self.path,
            )

        # Make these new models the output of this module
        self.output_models = renamed_models
        # Export outputs
        self.export_io_models()
