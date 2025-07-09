This guide explains how to run HADDOCK3 jobs either interactively or via SLURM batch submission (MPI HPC ready)
---

##  Quick Start

###  1. Run Interactively with `srun`

To  run HADDOCK3 in an interactive SLURM session:

#example code

srun --partition=compute \
     --nodes=1 \
     --ntasks-per-node=8 \
     --chdir=/path/to/haddock3/examples/docking-protein-protein \
     apptainer exec \
     --bind /path/to/host:/path/to/host \
     /path/to/haddock3_cpu-mpi.sif \
     haddock3 docking-protein-protein-mpi.cfg

 *Note: Change the `.cfg` file path and directory bindings to match your project location.*

---

###  2. Run as a Batch Job with `sbatch`

Submit your job to SLURM using the sample script in the `scripts/` folder:

```bash
scripts/ sbatch run_haddock3_slurm.sh
```

 *Customize the script and config paths as needed.*

---




##  Requirements

- SLURM-enabled HPC environment
- Docker, Apptainer or Singularity installed
- Image: `haddock-mpi-gpt.sif`
- `.cfg` configuration file for HADDOCK3

