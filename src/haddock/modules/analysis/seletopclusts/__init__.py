"""HADDOCK3 module to select a top cluster/model."""
import math
from pathlib import Path

from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """Haddock Module for 'seletopclusts'."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if module is installed."""
        return

    def _run(self):
        """Execute the module's protocol."""
        if self.params["top_models"] <= 0:
            _msg = "top_models must be either > 0 or nan."
            self.finish_with_error(_msg)

        if not isinstance(self.params["top_cluster"], list):
            _msg = "top_cluster must be a list, it can be an empty one."
            self.finish_with_error(_msg)

        models_to_select = self.previous_io.retrieve_models()

        # how many models should we output?
        self.output_models = []
        if not self.params["top_cluster"]:
            target_rankings = list(set([p.clt_rank for p in models_to_select]))
            target_rankins_str = ",".join(map(str, target_rankings))
            self.log(f"Selecting all clusters: {target_rankins_str}")
        else:
            target_rankings = list(set(self.params["top_cluster"]))
            target_rankins_str = ",".join(map(str, target_rankings))
            self.log(f"Selecting clusters: {target_rankins_str}")

        for target_ranking in target_rankings:
            if math.isnan(self.params["top_models"]):
                for pdb in models_to_select:
                    if pdb.clt_rank == target_ranking:
                        self.output_models.append(pdb)
            else:
                for model_ranking in range(1, self.params["top_models"] + 1):
                    for pdb in models_to_select:
                        if (
                                pdb.clt_rank == target_ranking
                                and pdb.clt_model_rank == model_ranking
                                ):
                            self.log(
                                f"Selecting {pdb.file_name} "
                                f"Cluster {target_ranking} "
                                f"Model {model_ranking}"
                                )
                            self.output_models.append(pdb)

        self.export_ouput_models()
