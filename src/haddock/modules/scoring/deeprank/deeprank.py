"""Wrapper around deeprank-gnn-esm for use as a HADDOCK3 scoring backend."""

import csv
import os
import sys
import tempfile
from pathlib import Path

from pdbtools import pdb_mkensemble


def deeprank_is_available() -> bool:
    """Check whether deeprank-gnn-esm's requirements are met."""
    try:
        import sqlite3  # noqa: F401
    except ImportError as err:
        raise ImportError(
            "The `deeprank` module requires a python interpreter built with "
            "sqlite3 support, which is missing from your current "
            "environment."
        ) from err

    try:
        import deeprank_gnn  # type: ignore  # noqa: F401
    except ImportError as err:
        raise ImportError(
            "The `deeprank` module requires the `deeprank_gnn` package, "
            "which is not installed in your current environment."
        ) from err

    return True


class DeeprankWrapper:
    """Run deeprank-gnn-esm on a set of models and collect predicted scores."""

    ENSEMBLE_NAME = "deeprank_ensemble.pdb"

    def __init__(self, models, ncores, chain_i, chain_j):
        self.models = models
        self.chain_i = chain_i
        self.chain_j = chain_j
        self.ncores = ncores

    def _make_ensemble(self, workspace: Path) -> Path:
        """Combine all input models into a single multi-model PDB.

        deeprank-gnn-esm only parallelizes (embedding batching, graph
        generation with `nproc`, and `NeuralNet(num_workers=...)`) across the
        models it finds inside one multi-MODEL PDB file.
        """
        ensemble_path = Path(workspace, self.ENSEMBLE_NAME)
        lines = pdb_mkensemble.run([str(model) for model in self.models])
        with open(ensemble_path, "w") as fh:
            fh.writelines(lines)
        return ensemble_path

    def run(self) -> dict[str, float]:
        """Score all models and return the predicted scores per model.

        deeprank-gnn-esm writes its ensemble/graph/embedding intermediates
        next to the input file with no way to redirect them, and none of
        that is useful once the scores are parsed out. So the whole run
        (ensemble build, scoring, csv parsing) happens inside a temporary
        directory that is discarded on exit, and only the scores dict is
        returned to the caller.
        """
        # This import needs to be exactly here
        from deeprank_gnn.predict import main as deeprank_main

        with tempfile.TemporaryDirectory() as workspace_str:
            workspace = Path(workspace_str)
            ensemble_path = self._make_ensemble(workspace)

            # NOTE: Since we are using the `main` function that takes `sys.argv`
            #  we need a hacky solution to override. Here we can simply re-write
            #  it and pass the arguments we need
            original_argv = sys.argv
            original_cwd = os.getcwd()
            sys.argv = [
                "deeprank",
                str(ensemble_path),
                self.chain_i,
                self.chain_j,
                str(self.ncores),
            ]

            try:
                # NOTE: deeprank will write its output to the path its being executed, there
                #  is no way to define where the output will be saved, so here we need to move
                #  into the workspace to trigger the function
                os.chdir(workspace)
                deeprank_main()
            finally:
                # NOTE: !!! VERY IMPORTANT !!!
                #  Since we moved directories and overrode the `sys.argv` we NEED to have this
                #  `finally` here - it means this branch of the code will always be executed.
                #  With this we can hopely guarantee we go back to where we should be before
                #  the execution moves on
                sys.argv = original_argv
                os.chdir(original_cwd)

            return self._retrieve_scores(workspace)

    def _retrieve_scores(self, workspace: Path) -> dict[str, float]:
        """Parse the output from deeprank and return the scores per input model.

        deeprank names each split model after its position in the ensemble
        file (`{ensemble_stem}_model_{index}`, 0-indexed in input order), so
        that naming is used to map predictions back onto `self.models`.
        """
        ensemble_stem = Path(self.ENSEMBLE_NAME).stem
        csv_path = (
            workspace
            / f"{ensemble_stem}-gnn_esm_pred_{self.chain_i}_{self.chain_j}"
            / "GNN_esm_prediction.csv"
        )
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            fnat_by_pdb_id = {
                row["pdb_id"]: float(row["predicted_fnat"]) for row in reader
            }

        scores = {}
        for index, model in enumerate(self.models):
            pdb_id = f"{ensemble_stem}_model_{index}"
            scores[str(model)] = fnat_by_pdb_id[pdb_id]
        return scores
