"""Wrapper around deeprank-gnn-esm for use as a HADDOCK3 scoring backend."""

import csv
import os
import shutil
import sys
import tempfile
from pathlib import Path


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

    INTERFACES_CSV = "GNN_esm_prediction.csv"
    SUMMARY_CSV = "GNN_esm_prediction_summary.csv"

    def __init__(self, models, ncores, score_type="mean", chain_i="A", chain_j="B"):
        self.models = models
        self.ncores = ncores
        self.score_type = score_type
        self.chain_i = chain_i
        self.chain_j = chain_j

    def run(self, output_dir: Path) -> dict[str, float]:
        """Score all models and return the predicted fnat per model.

        deeprank-gnn-esm predicts the fnat of every chain pair (interface) of
        each model and writes the per-interface and per-model csv files in the
        directory it is executed in. So it runs inside a temporary directory
        that is discarded on exit, and only these two csv files are copied to
        `output_dir`.
        """
        # This import needs to be exactly here
        from deeprank_gnn.main import main as deeprank_main

        # absolute paths, as deeprank is executed in another directory
        model_paths = [str(Path(model).resolve()) for model in self.models]

        with tempfile.TemporaryDirectory() as workspace_str:
            workspace = Path(workspace_str)

            # NOTE: Since we are using the `main` function that takes `sys.argv`
            #  we need a hacky solution to override. Here we can simply re-write
            #  it and pass the arguments we need
            original_argv = sys.argv
            original_cwd = os.getcwd()
            sys.argv = ["deeprank", *model_paths, "--num_cores", str(self.ncores)]

            try:
                # NOTE: deeprank writes its csv files in the path its being executed,
                #  so here we need to move into the workspace to trigger the function
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

            for fname in (self.INTERFACES_CSV, self.SUMMARY_CSV):
                shutil.copy(workspace / fname, Path(output_dir, fname))

            return self._retrieve_scores(workspace)

    def _retrieve_scores(self, workspace: Path) -> dict[str, float]:
        """Parse the output from deeprank and return the score per input model.

        deeprank names each model after its input file name without extension.
        The score is the mean or minimum predicted fnat over the interfaces of
        a model (from the summary csv), or the predicted fnat of the
        `chain_i`-`chain_j` interface (from the interfaces csv). An interface
        without contacts, or a model without any, has no prediction (`NA`) and
        gets a fnat of 0.
        """
        if self.score_type == "pair":
            with open(Path(workspace, self.INTERFACES_CSV)) as f:
                fnat_by_pdb_id = {
                    row["pdb_id"]: row["predicted_fnat"]
                    for row in csv.DictReader(f)
                    # deeprank sorts the chain IDs of each pair
                    if (row["chain_i"], row["chain_j"])
                    == tuple(sorted((self.chain_i, self.chain_j)))
                }
        else:
            with open(Path(workspace, self.SUMMARY_CSV)) as f:
                fnat_by_pdb_id = {
                    row["pdb_id"]: row[f"{self.score_type}_fnat"]
                    for row in csv.DictReader(f)
                }

        scores = {}
        for model in self.models:
            pdb_id = Path(model).stem
            if pdb_id not in fnat_by_pdb_id:
                raise ValueError(
                    f"Chains {self.chain_i} and {self.chain_j} not found "
                    f"in model {model}."
                )
            fnat = fnat_by_pdb_id[pdb_id]
            scores[str(model)] = 0.0 if fnat == "NA" else float(fnat)
        return scores
