"""HADDOCK3 module to select a top cluster/model."""
from pathlib import Path

from haddock.libs.libontology import ModuleIO
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yml")


class HaddockModule(BaseHaddockModule):
    """Haddock Module for 'seletopclusts'."""

    name = RECIPE_PATH.name

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
        """Confirm if module is installed."""
        return

    def _run(self):
        """Execute the module's protocol."""
        models_to_select = self.previous_io.retrieve_models()

        # how many models should we output?
        models = []
        for target_ranking in self.params['top_cluster']:
            if self.params['top_models'] == 'all':
                for pdb in models_to_select:
                    if pdb.clt_rank == target_ranking:
                        models.append(pdb)
            else:
                for model_ranking in range(1, self.params['top_models'] + 1):
                    for pdb in models_to_select:
                        if (pdb.clt_rank == target_ranking
                                and pdb.clt_model_rank == model_ranking):
                            self.log(
                                f"Selecting {pdb.file_name} "
                                f"Cluster {target_ranking} "
                                f"Model {model_ranking}")
                            models.append(pdb)

        # Save module information
        io = ModuleIO()
        io.add(models, "o")
        io.save()
