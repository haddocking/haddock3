import logging
import shutil
from os import linesep
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.ontology import Format, ModuleIO, PDBFile
from haddock.pdbutil import PDBFactory


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

        # Check if multiple models are provided
        if len(models_to_score) > 1:
            self.finish_with_error("Only one model allowed in LightDock sampling module")

        model = models_to_score[0]
        # Check if chain IDs are present
        segids, chains = PDBFactory.identify_chainseg(Path(model.path) / model.file_name)
        if set(segids) != set(chains):
            logger.info("No chain IDs found, using segid information")
            PDBFactory.swap_segid_chain(Path(model.path) / model.file_name,
                self.path / model.file_name)
        else:
            # Copy original model to this working path
            shutil.copyfile(Path(model.path) / model.file_name, self.path / model.file_name)

        model_with_chains = self.path / model.file_name
        # Split by chain
        new_models = PDBFactory.split_by_chain(model_with_chains)
        print(new_models)

        # Dummy expected
        expected = models_to_score

        # Save module information
        io = ModuleIO()
        for model in models_to_score:
            io.add(model)
        for p in expected:
            io.add(p, "o")
        io.save(self.path)
