#!/bin/bash

#SBATCH --job-name=HADDOCK3-docking
#SBATCH --output=HADDOCK3_%j.out
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=120
#SBATCH --mem=256GB
#SBATCH --partition=compute

echo "Starting HADDOCK3 Docking Job"
echo "SLURM_JOBID          = $SLURM_JOBID"
echo "SLURM_JOB_NODELIST   = $SLURM_JOB_NODELIST"
echo "SLURM_NNODES         = $SLURM_NNODES"
echo "SLURMTMPDIR          = $SLURMTMPDIR"
echo "Date                 = $(date)"
echo "Hostname             = $(hostname -s)"
echo "Working Directory    = $(pwd)"
echo "Submit Directory     = $SLURM_SUBMIT_DIR"

# Load necessary environment (optional depending on your setup)
source /lustre/oneApi/setvars.sh
export OMP_NUM_THREADS=1

# Run HADDOCK3 via Apptainer
cd "${WORK_DIR}"

# Ensure the working directory is set correctly
apptainer exec --bind /path/to/host/data:/path/to/container/data \
  /path/to/haddock3_image.sif \
  haddock3 your-docking-config.cfg


echo "HADDOCK3 Job Complete"
echo "Completed at: $(date)"
