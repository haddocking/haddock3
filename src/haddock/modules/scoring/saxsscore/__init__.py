"""SAXS scoring module."""
from pathlib import Path
import subprocess
from os import linesep
import numpy as np

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
        ns = self.params["ns"]
        cst = self.params["cst"]
        w_haddock = self.params["w_haddock"]
        w_saxs = self.params["w_saxs"]

        for model in models_to_score:

            # Calculate chi value
            cmd = f"crysol {model.rel_path} {saxs_data} -lm {lm} -ns {ns}"
            if cst:
                cmd += " -cst"
            subprocess.call(cmd, shell=True)

            # Read chi value
            fit_file = model.file_name.replace(".pdb", "00.fit")
            with open(fit_file, "r") as ff:
                fit_header = ff.readline()
                model.chi = float(fit_header.partition("Chi^2:")[2].partition(r"\s")[0])  # noqa: E501

            # Calculate HADDOCKsaxs score
            if np.isnan(model.score):
                self.log("No HADDOCK score from previous module")
                self.log("Using chi value as score")
                haddock_score = model.chi
            else:
                self.log("Calculating HADDOCKsaxs score")
                haddock_score = (model.score * w_haddock) + (model.chi * w_saxs)

            model.score = haddock_score

        self.output_models = models_to_score

        # Write scores to output file
        output_fname = "saxsscore.tsv"
        self.log(f"Saving output to {output_fname}")

        text_generator = (
            f"{pdb.file_name}\t{pdb.ori_name}\t{pdb.md5}\t{pdb.chi}\t{pdb.score}"  # noqa: E501
            for pdb in self.output_models)

        with open(output_fname, "w") as fh:
            fh.write("\t".join(
                ("structure", "original_name", "md5", "chi", "score")
                ) + linesep)
            fh.write(linesep.join(text_generator))

        # Send models (unchanged but with HADDOCKsaxs score) to the next step
        self.export_output_models()
