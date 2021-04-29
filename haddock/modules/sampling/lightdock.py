import logging
from os import linesep
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.ontology import Format, ModuleIO, PDBFile


logger = logging.getLogger(__name__)


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path):
        recipe_path = Path(__file__).resolve().parent.absolute()
        defaults = recipe_path / "sampling.toml"
        super().__init__(order, path, defaults=defaults)

    def run(self, module_information):
        logger.info("Running [sampling-lightdock] module")

        # Apply module information to defaults
        self.patch_defaults(module_information)

        # Get the models generated in previous step
        models_to_score = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        # Dummy expected
        expected = models_to_score

        # Save module information
        io = ModuleIO()
        for model in models_to_score:
            io.add(model)
        for p in expected:
            io.add(p, "o")
        io.save(self.path)
