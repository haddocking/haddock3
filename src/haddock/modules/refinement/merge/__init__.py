"""HADDOCK3 module to select a top cluster/model."""
from itertools import chain, product
from pathlib import Path

from haddock import log
from haddock.libs import libpdb
from haddock.libs.libontology import Format, ModuleIO, PDBFile
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to merge topologies into complexes."""

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if module is installed."""
        return

    def run(self, **params):
        """Execute module."""
        log.info("Running [merge] module")

        super().run(params)

        # Get the models generated in previous step
        #  The models from topology come in a dictionary
        input_dic = {}
        for i, model_dic in enumerate(self.previous_io.output):
            input_dic[i] = []
            for key in model_dic:
                input_dic[i].append(model_dic[key])

        # Generate all combinations between the ensembles
        models_to_be_merged = [
            values for values in product(*input_dic.values())
            ]
        merged_models = []
        for i, model_list in enumerate(models_to_be_merged, start=1):
            pdb_fname = Path(self.path, f"merged_{i}.{Format.PDB}")
            pdb_list = [e.full_name for e in model_list]

            libpdb.merge(pdb_list, pdb_fname)

            topology = list(chain(*[e.topology for e in model_list]))

            pdb = PDBFile(pdb_fname, topology, path=self.path)

            merged_models.append(pdb)

        io = ModuleIO()
        io.add(merged_models, "o")
        io.save(self.path)
