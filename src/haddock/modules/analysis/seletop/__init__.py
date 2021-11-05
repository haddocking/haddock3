"""HADDOCK3 module to select a top cluster/model"""
import logging
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.libs.libontology import Format, ModuleIO

logger = logging.getLogger(__name__)

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):

    def __init__(
            self,
            order,
            path,
            *ignore,
            init_params=DEFAULT_CONFIG,
            **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        return

    def run(self, **params):
        logger.info("Running [seletop] module")

        super().run(params)

        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            self.finish_with_error('This module cannot come after one'
                                   ' that produced an iterable')

        models_to_select = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        # sort the models based on their score
        models_to_select.sort(key=lambda x: x.score)

        if len(models_to_select) < self.params['select']:
            logger.warning('Number of models to be selected is larger'
                           ' than generated models, selecting ALL')

        # select the models based on the parameter
        selected_models = models_to_select[:self.params['select']]

        io = ModuleIO()
        io.add(selected_models, "o")
        io.save(self.path)
