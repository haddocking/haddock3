# HADDOCK3 HPC Usage and Configuration Guide

This guide explains how to run HADDOCK3 jobs either interactively or via SLURM batch submission (MPI HPC ready)
---

## Quick Start

### 1. Run interactively with `srun`

To run HADDOCK3 in an interactive SLURM session:

#example code

srun --partition=compute \
     --nodes=1 \
     --ntasks-per-node=8 \
     --chdir=/path/to/haddock3/examples/docking-protein-protein \
     /path/to/haddock3_cpu-mpi.sif \
     haddock3 docking-protein-protein-mpi.cfg

<pre>
<strong>Note:</strong> Adjust paths— change the `.cfg` file path and directory to match your file location.
</pre>



<pre>
<strong>CPU allocation and container behavior:</strong> No rebuild is required for the container — it will automatically use the number of CPUs assigned by SLURM.
In the above example, `--nodes=1` allocates one compute node, and `--ntasks-per-node=8` allocates eight parallel tasks. The exact number of nodes and CPU cores you can request depends on your HPC  configuration — please refer to your cluster documentation for resource allocation limit.
</pre>

---

### 2. Run as a batch job with `sbatch`

Submit your job to SLURM using the sample script in the `scripts/` folder:

```bash
scripts/sbatch run_haddock3_slurm.sh
```

*Customize the script and configuration paths as needed.*

---

### 3. MPI run mode

When using multiple CPUs or running across multiple nodes, **HADDOCK3 must be executed in MPI mode**.
To ensure proper parallel execution, set the following parameter in your configuration (`.cfg`) file, and change mode and the number of ncores:

```ini
mode = "mpi"
ncores = 8
```

This enables HADDOCK3 to distribute workloads efficiently across all allocated CPUs.

---

## 4. Requirements

- SLURM-enabled HPC environment
- Docker, Apptainer installed
- Image: `haddock3_cpu-mpi.sif`
- `.toml` configuration file for HADDOCK3
- Root (sudo) privileges are required to build Apptainer images locally
