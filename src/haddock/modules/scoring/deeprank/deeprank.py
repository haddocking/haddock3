"""Wrapper around deeprank-gnn-esm for use as a HADDOCK3 scoring backend."""

import csv
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

    def __init__(self, models, ncores, score_type, chain_i, chain_j):
        self.models = models
        self.ncores = ncores
        self.score_type = score_type
        self.chain_i = chain_i
        self.chain_j = chain_j

    def run(self) -> dict[str, float]:
        """Score all models and return the predicted fnat per model.

        deeprank-gnn-esm predicts the fnat of every chain pair (interface) of
        each model. It keeps its intermediate files in a temporary directory
        that it removes, and writes only the per-interface and per-model csv
        files in the current directory, from where the scores are read.
        """
        # This import needs to be exactly here
        from deeprank_gnn.main import run_predictions

        run_predictions(self.models, num_cores=self.ncores)

        return self._retrieve_scores(Path.cwd())

    def _retrieve_scores(self, directory: Path) -> dict[str, float]:
        """Parse the output from deeprank and return the score per input model.

        deeprank names each model after its input file name without extension.
        The score is the mean or minimum predicted fnat over the interfaces of
        a model (from the summary csv), or the predicted fnat of the
        `chain_i`-`chain_j` interface (from the interfaces csv). An interface
        without contacts, or a model without any, has no prediction (`NA`) and
        gets a fnat of 0.
        """
        if self.score_type == "pair":
            # deeprank sorts the chain IDs of each pair
            pair = tuple(sorted((self.chain_i, self.chain_j)))
            csv_name, column = self.INTERFACES_CSV, "predicted_fnat"
        else:
            pair = None
            csv_name, column = self.SUMMARY_CSV, f"{self.score_type}_fnat"

        with open(Path(directory, csv_name)) as f:
            fnat_by_pdb_id = {
                row["pdb_id"]: row[column]
                for row in csv.DictReader(f)
                if pair is None or (row["chain_i"], row["chain_j"]) == pair
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
