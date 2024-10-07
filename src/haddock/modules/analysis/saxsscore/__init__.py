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
from pathlib import Path

from haddock.libs.libutil import check_subprocess
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.saxsscore.saxs import (
    run_crysol,
    read_chi2,
    calculate_haddocksaxs_score,
    output_saxsscore,
)


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
        check_subprocess("crysol -h")

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

            # Run CRYSOL
            run_crysol(
                atsas_path=self.params["atsas_path"],
                input_f=model.rel_path,
                saxs_data=saxs_data,
                lm=lm,
                ns=ns,
                cst=cst,
            )

            # Read Chi^2 value
            fit_file = (Path(model.file_name).stem 
                        + "_" + Path(saxs_data).stem + ".fit")
            model.chi2 = read_chi2(fit_file=fit_file)

            # Calculate HADDOCKsaxs score
            if math.isnan(model.score):
                haddocksaxs_score = model.chi2
            else:
                haddocksaxs_score = calculate_haddocksaxs_score(
                    score=model.score,
                    chi2=model.chi2,
                    w_haddock=w_haddock,
                    w_saxs=w_saxs,
                )

            model.score = haddocksaxs_score

        # Write output
        self.output_models = models_to_score
        output_fname="saxsscore.tsv"
        self.log(f"Saving output to {output_fname}")

        output_saxsscore(
            output_models=self.output_models,
            output_fname=output_fname,
        )

        # Send models (unchanged) to the next step
        self.export_io_models()
