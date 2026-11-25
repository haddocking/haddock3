# HADDOCK3 examples

Each folder contains example configuration files for one scenario. Files
ending in `-test.cfg` are short test runs; `-full.cfg` files are production-like
runs. Run an example from inside its folder:

```bash
cd docking-protein-protein
haddock3 docking-protein-protein-test.cfg
```

# Running jobs with MPI

See our documentation page at `docs/pages/tutorials/mpi.md`.

# Regression testing against the main branch

`run_examples_regression.py` runs every `*/*-test.cfg` example, checks that each
run completed without errors, and compares the results with a reference built
from `origin/main`.

```bash
python run_examples_regression.py --list    # selected and skipped tests
python run_examples_regression.py init      # build the reference from main
python run_examples_regression.py run       # run the current code and compare
python run_examples_regression.py compare   # compare saved results only
```

- `init` checks out `origin/main` in a separate git worktree
  (`../haddock3-main-ref` relative to the repository, see `--ref-worktree`)
  and runs the examples there; your checkout is not touched.
- Examples using MPI or batch submission, and follow-up configurations that
  need an existing run (`--extend-run`, `--restart`), are skipped.
- A run is considered successful when haddock3 exits without error, the run log
  ends with `Finished at` and contains no `ERROR`/`CRITICAL` lines, and every
  module has its step folder. The number of output models and the files written
  by each step are compared with the reference.
- For every caprieval step, the CAPRI statistics (number of acceptable, medium
  and high quality models, rank of the first acceptable model, rank and quality
  of the best model) are compared. These are the statistics of the
  `extract-capri-stats-dockq.sh` script (and `extract-capri-stats-ligand.sh`
  for the protein-ligand examples).
- Use `--only PATTERN` (repeatable) to run a subset, e.g. `--only scoring/`.
  Merging with saved results from another commit is refused unless
  `--allow-mixed` is given.
  `-h` lists the tolerances and other options.
- Results, logs and the report (`report.md`) are written to
  `examples/regression/`.

Running all examples takes several hours, depending on your hardware.

**Warning:** existing run directories of the examples are deleted before each
run.

# Extracting CAPRI statistics from a run

The `extract-capri-stats-*.sh` scripts summarise the CAPRI tables written by
each `caprieval` step of a run. Give them the run directory as argument:

```bash
./extract-capri-stats-dockq.sh docking-protein-protein/run1-test
```

```
==============================================
== docking-protein-protein/run1-test/2_caprieval/capri_ss.tsv
==============================================
Total number of acceptable or better models:  11  out of  20
Total number of medium or better models:      5  out of  20
Total number of high quality models:          0  out of  20

First acceptable model - rank:  1  i-RMSD:  2.207  Fnat:  0.389  DockQ:  0.473
First medium model     - rank:  2  i-RMSD:  1.372  Fnat:  0.417  DockQ:  0.621
Best model             - rank:  4  i-RMSD:  0.992  Fnat:  0.667  DockQ:  0.762
...
```

For each caprieval step they report the number of acceptable or better, medium
or better and high quality models (or clusters), the first acceptable and first
medium model in the ranking, and the best model. Lower i-RMSD is better, higher
DockQ is better.

| Script | Table | Quality criteria (acceptable / medium / high) | Best model |
|---|---|---|---|
| `extract-capri-stats-dockq.sh` | models (`capri_ss.tsv`) | DockQ ≥ 0.23 / 0.49 / 0.80 | lowest i-RMSD |
| `extract-capri-stats-irmsd.sh` | models (`capri_ss.tsv`) | i-RMSD ≤ 4 / 2 / 1 Å | lowest i-RMSD |
| `extract-capri-stats-ligand.sh` | models (`capri_ss.tsv`) | il-RMSD ≤ 2.5 / 1.0 / 0.5 Å | lowest il-RMSD |
| `extract-capri-stats-dockq-clt.sh` | clusters (`capri_clt.tsv`) | mean DockQ ≥ 0.23 / 0.49 / 0.80 | lowest mean i-RMSD |
| `extract-capri-stats-irmsd-clt.sh` | clusters (`capri_clt.tsv`) | mean i-RMSD ≤ 4 / 2 / 1 Å | lowest mean i-RMSD |

The statistics are only meaningful when a reference structure is given to
`caprieval` (`reference_fname`); without it, the best-scoring model is used as
reference. Use the ligand script for small
molecule docking, where the il-RMSD (ligand RMSD after fitting on the
interface) is the relevant measure.
