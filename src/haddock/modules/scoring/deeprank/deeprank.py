import csv
import os
import sys
from pathlib import Path


def deeprank_is_available() -> bool:
    try:
        import deeprank_gnn  # type: ignore
    except ImportError:
        raise
    return True


class DeeprankWraper:
    def __init__(self, models, ncores, chain_i, chain_j, path):
        self.models = models
        self.chain_i = chain_i
        self.chain_j = chain_j
        self.ncores = ncores
        self.path = path

    def run(self):
        """Run method for the wrapper, it will call deeprank as if we were using the `main` function."""

        # This import needs to be exactly here
        from deeprank_gnn.predict import main as deeprank_main

        for model in self.models:
            # NOTE: Since we are using the `main` function that takes `sys.argv`
            #  we need a hacky solution to override. Here we can simply re-write
            #  it and pass the arguments we need
            original_argv = sys.argv
            original_cwd = os.getcwd()
            sys.argv = [
                "deeprank",
                str(model),
                self.chain_i,
                self.chain_j,
                str(self.ncores),
            ]

            try:
                # NOTE: deeprank will write its output to the path its being executed, there
                #  is no way to define where the output will be saved, so here we need to move
                #  into the `self.path` to trigger the function
                os.chdir(self.path)
                deeprank_main()
            finally:
                # NOTE: !!! VERY IMPORTANT !!!
                #  Since we moved directories and overrode the `sys.argv` we NEED to have this
                #  `finally` here - it means this branch of the code will always be executed.
                #  With this we can hopely guarantee we go back to where we should be before
                #  the execution moves on
                sys.argv = original_argv
                os.chdir(original_cwd)

    def retrieve_scores(self) -> dict[str, float]:
        """Parse the output from deeprank and return the scores."""
        scores = {}
        for model in self.models:
            csv_path = (
                Path(self.path)
                / f"{Path(model).stem}-gnn_esm_pred_{self.chain_i}_{self.chain_j}"
                / "GNN_esm_prediction.csv"
            )
            with open(csv_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    scores[str(model)] = float(row["predicted_fnat"])
        return scores
