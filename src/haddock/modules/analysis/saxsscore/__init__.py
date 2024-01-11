"""
SAXS scoring module
===================

Requires ATSAS 2.6.0 or newer.

This module adds a SAXS term to the HADDOCK score.
The SAXS term is based on the fit (Chi^2) of a model to corresponding
scattering data, as calculated by the ATSAS program CRYSOL.

By default the weights of the HADDOCK score (calculated by the preceding module)
and the Chi value are 1 and 50, respectively.
If the preceding module did not provide HADDOCK scores,
the Chi^2 values are returned as is (not multiplied by 50).

The resulting scores are found in ``saxsscore.tsv``.
This module does not alter the models themselves
and passes them on unchanged to the subsequent module.

For more information on SAXS-based scoring in HADDOCK, see

Karaca E, Bonvin AMJJ. 2013. On the Usefulness of Ion-Mobility Mass Spectrometry
and SAXS Data in Scoring Docking Decoys. *Acta Crystallographica* D69:683-694.
"""
import math
import subprocess
from os import linesep
from pathlib import Path

from haddock.libs.libutil import check_subprocess
from haddock.modules import BaseHaddockModule


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

            # Calculate Chi^2 value
            cmd = f"crysol {model.rel_path} {saxs_data} -lm {lm} -ns {ns}"
            if cst:
                cmd += " -cst"
            subprocess.call(cmd, shell=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            # Read Chi^2 value
            fit_file = model.file_name.replace(".pdb", "00.fit")
            with open(fit_file, "r") as ff:
                fit_header = ff.readline()
                model.chi2 = float(fit_header.partition("Chi^2:")[2].partition(r"\s")[0])  # noqa: E501

            # Calculate HADDOCKsaxs score
            if math.isnan(model.score):
                self.log(f"Using Chi^2 value as score for {model.file_name}")
                haddock_score = model.chi2
            else:
                chi = math.sqrt(model.chi2)
                haddock_score = (model.score * w_haddock) + (chi * w_saxs)

            model.score = haddock_score

            fit_log = model.file_name.replace(".pdb", "00.log")
            cmd = f"rm {fit_file} {fit_log}"
            subprocess.call(cmd, shell=True)

        self.output_models = models_to_score

        # Write scores to output file
        output_fname = "saxsscore.tsv"
        self.log(f"Saving output to {output_fname}")

        text_generator = (
            f"{pdb.file_name}\t{pdb.ori_name}\t{pdb.md5}\t{pdb.chi2}\t{pdb.score}"  # noqa: E501
            for pdb in self.output_models)

        with open(output_fname, "w") as fh:
            fh.write("\t".join(
                ("structure", "original_name", "md5", "Chi^2", "score")
                ) + linesep)
            fh.write(linesep.join(text_generator))

        # Send models (unchanged) to the next step
        self.export_output_models()
