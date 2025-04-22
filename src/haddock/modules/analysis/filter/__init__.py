"""Filter models based on cutoff value.

This module filters models from the input models based on a cutoff value.
By default, the selection is based on the HADDOCK score of the models.

The number of models to be selected is unknown and by the parameter `select`.
In the standard HADDOCK protocol, this number is 200, which can be increased
if more models should be refined.
"""
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.libs.libontology import Format, PDBFile
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to select top cluster/model."""

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
        """Execute module."""
        # Make sure we have access to complexes
        if type(self.previous_io) == iter:
            self.finish_with_error(
                "[filter] This module cannot come after one"
                " that produced an iterable."
                )
        # Get the models generated in previous step
        models: list[PDBFile] = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.PDB
            ]
        
        # Get the filter by parameter
        filter_by = "score"
        threshold = self.params["cutoff"]

        # Make sure we can access this attribute on models
        without_score_values = {"FAKE", None}
        models_with_attributes: list[PDBFile] = [
            m for m in models
            if getattr(m, filter_by, "FAKE") not in without_score_values
            ]
        
        # Check how many of them are available
        ratio_models_with_attr = len(models_with_attributes) / len(models)
        self.log(
            f"{100 * (1 - ratio_models_with_attr):6.2f} % "
            "of the input models have accessible scores."
            )
        if len(models_with_attributes) == 0:
            self.finish_with_error(
                "Input models do not have scores. "
                "Please consider running a scoring module before!"
                )

        # Process to the actual filtering step
        filtered_models: list[PDBFile] = [
            m for m in models_with_attributes
            if getattr(m, filter_by) <= threshold
            ]

        # Final evaluation of the outcome of the filtering
        percent_filtered = (1 - (len(filtered_models) / len(models))) * 100
        if len(filtered_models) == 0:
            self.finish_with_error(
                f"With the currently set 'cutoff' value of {threshold}, "
                "ALL models were filtered out."
                )
        else:
            self.log(
                f"With currently set 'cutoff' value of {threshold}, "
                f"{percent_filtered:6.2f}% of the models were filtered out."
                )

        # select the models based on the parameter
        self.output_models = filtered_models
        self.export_io_models()
