"""Select a number of models.

This module selects a number of models from the input models. By default, the
selection is based on the HADDOCK score of the models.

The number of models to be selected is defined by the parameter `select`.
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
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            self.finish_with_error('[seletop] This module cannot come after one'
                                   ' that produced an iterable')

        models_to_select: list[PDBFile] = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.PDB
            ]

        # sort the models based on their score
        models_to_select.sort(key=lambda x: x.score)

        if len(models_to_select) < self.params['select']:
            self.log((
                'Number of models to be selected is larger'
                ' than generated models, selecting ALL'),
                level='warning',
                )

        # select the models based on the parameter
        self.output_models = models_to_select[:self.params['select']]
        self.export_io_models()
