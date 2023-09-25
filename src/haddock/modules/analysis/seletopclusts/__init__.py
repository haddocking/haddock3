"""Select a top cluster module."""
import math
import os
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

        if not isinstance(self.params["top_cluster"], int):
            _msg = "top_cluster must be an integer."
            self.finish_with_error(_msg)

        models_to_select = self.previous_io.retrieve_models()

        # how many models should we output?
        self.output_models = []
        target_rankings = list(set([p.clt_rank for p in models_to_select]))
        if self.params["top_cluster"] >= len(target_rankings):
            # select all clusters
            target_rankins_str = ",".join(map(str, target_rankings))
            self.log(f"Selecting all clusters: {target_rankins_str}")
        else:
            # select top_cluster clusters
            target_rankings = list(range(1, self.params["top_cluster"] + 1))
            target_rankins_str = ",".join(map(str, target_rankings))
            self.log((f"Selecting top {self.params['top_cluster']} clusters: "
                      f"{target_rankins_str}"))

        # select the models to export
        models_to_export = []
        for target_ranking in target_rankings:
            if math.isnan(self.params["top_models"]):
                for pdb in models_to_select:
                    if pdb.clt_rank == target_ranking:
                        models_to_export.append(pdb)
            else:
                for model_ranking in range(1, self.params["top_models"] + 1):
                    for pdb in models_to_select:
                        if (
                                pdb.clt_rank == target_ranking
                                and pdb.clt_model_rank == model_ranking
                                ):
                            self.log(
                                f" {pdb.file_name} "
                                f"> cluster_{target_ranking}_"
                                f"model_{model_ranking}.pdb"
                                )
                            models_to_export.append(pdb)
        
        # dump the models to disk and change their attributes
        with open('seletopclusts.txt', 'w') as fh:
            fh.write("rel_path\tori_name\tcluster_name\tmd5" + os.linesep)
            for model in models_to_export:
                name = (
                    f"cluster_{model.clt_rank}_model"
                    f"_{model.clt_model_rank}.pdb")
                # writing name
                fh.write(
                    f"{model.rel_path}\t"
                    f"{model.ori_name}\t"
                    f"{name}\t"
                    f"{model.md5}" + os.linesep)
                # changing attributes
                name_path = Path(name)
                name_path.write_text(model.rel_path.read_text())
                model.ori_name = model.file_name
                model.file_name = name
                model.full_name = name
                model.rel_path = Path('..', Path(self.path).name, name)
                model.path = str(Path(self.path).resolve())

        self.output_models = models_to_export
        self.export_io_models()
