# Running HADDOCK3 with MPI

To run this example you must have `mpi4py` installed in the haddock3 python
environment, and `OpenMPI` in the host system.

```bash
$ pip install mpi4py
# or
$ conda install -c conda-forge mpi4py
```

**Do not run it as: `mpirun -np haddock3 run.cfg`**, this logic is handled
internally.

Edit the `.cfg` according to the `.job` parameters you have set

In the `.job` header:

```bash
#SBATCH --nodes=5
#SBATCH --tasks-per-node=96
```

In the `.cfg` params

```toml
ncores = 480
```

Then prepare the rest of the `.job` file according to your cluster, a SLURM
example is provided at
`haddock3/examples/docking-protein-protein/docking-protein-protein-mpi.job` and
below:

```bash
#!/bin/bash
#SBATCH --nodes=5
#SBATCH --tasks-per-node=96
#SBATCH -J haddock3mpi

# make sure anaconda is activated
source $HOME/software/miniconda3/bin/activate
conda activate haddock3

# go to the example directory
cd $HOME/repos/haddock3/examples/docking-protein-protein

# remove any old runs
rm -rf run-mpi

# execute
haddock3 docking-protein-protein-mpi-test.cfg
```

Then:

```bash
sbatch docking-protein-protein-mpi.job
```

Please report any errors as an [issue](https://github.com/haddocking/haddock3/issues/new/choose).
