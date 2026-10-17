# Running the example regression tests

*This documentation is for developers.*

The example regression tests check that HADDOCK3 workflows still run correctly
and give the same results as the `main` branch. They run all short test
examples (`examples/*/*-test.cfg`) on your branch and on `main`, and compare:

- whether each run completed without errors;
- the outputs written by each step (number of models, output files);
- the CAPRI statistics of every `caprieval` step: number of acceptable, medium
  and high quality models, rank of the first acceptable model, and rank and
  quality of the best model.

The full set takes about 3 hours per branch on a local machine, depending on
the hardware. It is therefore not part of the CI. Please run it before
submitting a pull request that changes the core workflow, the CNS scripts or the
force field. If you do not have the compute power to run these tests, ask us and
we will run them for you.

## 1. Build the reference from `main`

From the `examples/` folder of your development repository, with the haddock3
python environment activated:

```bash
cd examples/
python run_examples_regression.py init
```

This creates (or updates) a git worktree of `origin/main` next to your
repository (`../haddock3-main-ref`, see `--ref-worktree`), runs the examples
there and saves the results to `examples/regression/reference.json`. Your own
checkout is not touched. The reference only needs to be rebuilt when `main`
changes.

## 2. Run the examples on your branch and compare

```bash
python run_examples_regression.py run
```

This runs the same examples with your code, saves the results to
`examples/regression/current.json` and writes the comparison to
`examples/regression/report.md`. The report lists:

1. tests that fail on your branch but not on `main` (with the reason and the
   end of the log), and steps whose outputs differ;
2. changes in the CAPRI statistics;
3. the rank of the first acceptable model and the rank and quality of the best
   model for every caprieval step, with a better/worse verdict;
4. significant runtime changes;
5. tests failing on both branches (e.g. missing third-party packages);
6. unchanged and skipped tests.

The runs use fixed random seeds, so unchanged code gives identical results. If
there are differences, explain them in the pull request description unless
something went wrong.

## Useful options

- `--list` shows which examples will run, which are skipped and why, and which
  existing run directories will be deleted.
- `--only PATTERN` (repeatable) runs a subset, e.g.
  `--only docking-protein-DNA --only scoring/`. The results are merged with
  those already saved, so you can rerun a single example and compare again.
- `compare` rewrites the report from the saved results without running
  anything, e.g. with other tolerances (`--tol-dockq`, `--tol-irmsd`,
  `--tol-runtime`).
- `-b` stops at the first failing example; `--timeout` limits the time per
  example.

Examples using MPI or batch submission, and follow-up configurations that need
an existing run (`--extend-run`, `--restart`), are skipped.

**Warning:** existing run directories of the examples are deleted before each
run.

## Inspecting a single run

The `examples/extract-capri-stats-*.sh` scripts print the CAPRI statistics of
all caprieval steps of a run directory, see `examples/README.md`:

```bash
./extract-capri-stats-dockq.sh docking-protein-protein/run1-test
```

*Thanks, and enjoy developing and using HADDOCK3 :-)*
