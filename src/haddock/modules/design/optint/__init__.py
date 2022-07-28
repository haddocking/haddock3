"""interface optimization ."""
import importlib
import sys
from pathlib import Path

from haddock import log
from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to perform energy minimization scoring."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm module is installed."""
        dependency_found = importlib.util.find_spec("deap") is not None
        if not dependency_found:
            log.error(
                "Cannot execute `optint` module, "
                "install the dependency with `pip install deap`"
                )
            sys.exit(1)

    def _run(self):
        """Execute module."""
        from .optint import run

        try:
            input_pdb_l = self.previous_io.retrieve_models()
        except Exception as e:
            self.finish_with_error(e)

        output_model_list = run(input_pdb_l, self.params)

        self.output_models = []
        for element in output_model_list:
            pdb, psf, score = element
            # IMPORTANT: pass `file_name=Path.name` or CNS will fail
            pdb_object = PDBFile(
                Path(pdb).name,
                topology=TopologyFile(
                    Path(psf).name, path="."
                    ),
                path=".")
            pdb_object.score = score

            self.output_models.append(pdb_object)

        self.export_output_models()
