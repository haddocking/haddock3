"""SAXS scoring module."""
from pathlib import Path
import subprocess

from haddock.gear.haddockmodel import HaddockModel
from haddock.modules import BaseHaddockModule
from haddock.libs.libutil import check_subprocess


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to perform SAXS-based scoring."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm this module is installed."""
        check_subprocess('crysol -h')

    def _run(self):
        """Execute module."""

        models_to_score = self.previous_io.retrieve_models(individualize=True)

        # Get parameters from .yaml/.cfg
        saxs_data = self.params["saxs_fname"]
        lm = self.params["lm"]
        cst = self.params["cst"]
        ns = self.params["ns"]

        for model_num, model in enumerate(models_to_score, start=1):

            # Calculate chi values
            cmd = f"crysol {model.rel_path} {saxs_data} -lm {lm} {cst} -ns {ns}"
            subprocess.call(cmd, shell=True)

        # Send models (unchanged) to the next step
        self.output_models = models_to_score
        self.export_output_models()
