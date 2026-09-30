# GPU Acceleration in HADDOCK3

HADDOCK3 includes optional, zero-regression GPU acceleration designed to eliminate the heaviest computational bottlenecks in biomolecular docking, refinement, and ensemble analysis.

---

## 1. Quick Start

Install the optional GPU dependencies in your HADDOCK3 virtual environment:

```bash
pip install 'haddock3[gpu]'
```

Enable GPU acceleration in your workflow configuration file (`run.cfg`):

```yaml
# Global workflow parameters
use_gpu: true
gpu_platform: "auto"       # Options: "auto", "cuda", "mps", "opencl", "cpu"
gpu_devices: [0]           # List of GPU device IDs to allocate
```

---

## 2. Configuration Parameters

The following parameters are supported globally in `run.cfg` or individually within specific modules:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `use_gpu` | `boolean` | `false` | Enables GPU acceleration across supported modules with automatic fallback to CPU. |
| `gpu_platform` | `string` | `"auto"` | Target platform backend: `"auto"` (prioritizes CUDA > MPS > OpenCL > CPU), `"cuda"`, `"mps"` (Apple Silicon), `"opencl"`, or `"cpu"`. |
| `gpu_devices` | `list[int]` | `[]` | Explicit list of GPU device indices to assign among parallel worker processes (e.g. `[0, 1, 2, 3]`). |

---

## 3. Accelerated Modules Overview

### 3.1 Pairwise RMSD Matrix (`rmsdmatrix` and `ilrmsdmatrix`)

In standard execution, pairwise RMSD analysis across N structural models scales as O(N^2) and requires exporting coordinates to disk (`traj.xyz`) for C execution.

With `use_gpu: true`:
- Coordinates are loaded into memory and aligned via batched Kabsch Singular Value Decomposition (SVD) in double precision (`float64`).
- Centering, covariance cross-products, and reflection corrections are evaluated entirely in GPU VRAM.
- Achieves 15x to 50x speedups over CPU execution for ensembles of 1,000+ structures.

### 3.2 Fraction of Common Contacts Clustering (`clustfcc`)

The Fraction of Common Contacts (FCC) algorithm calculates the residue-residue contact overlap between all pairs of decoy complexes. In pure Python, evaluating N * (N - 1) / 2 set intersections is a major bottleneck for large ensembles.

With `use_gpu: true`:
- Contacts are encoded as a binary occurrence matrix A of shape (N, total_unique_contacts).
- Pairwise intersection counts are evaluated simultaneously using matrix multiplication: M = A * A^T.
- Accelerated using GPU Tensor Cores or compiled SciPy sparse routines, reducing clustering preparation from minutes to milliseconds.

### 3.3 Inter-Chain Contact Heatmaps (`contactmap`)

Evaluating all-atom pairwise distances for large multi-protein complexes can exceed system host RAM.

With `use_gpu: true`:
- Uses chunked GPU `torch.cdist` in single precision (`float32`), reducing memory usage by 50%.
- Work is processed in row chunks of 2,048 atoms, guaranteeing that VRAM limits are respected even for complexes exceeding 50,000 atoms.
- Bypasses the defensive host RAM abort guard.

### 3.4 Molecular Dynamics Refinement (`openmm`)

The OpenMM refinement module performs energy minimization and short molecular dynamics in implicit or explicit solvent.
- Automatically selects the `CUDA` or `OpenCL` platform when `use_gpu: true`.
- Sets `Precision: mixed` for optimal throughput on modern Tensor Core hardware.
- Distributes independent refinement trajectories across available GPUs using process-level affinity.

### 3.5 Scoring and Sampling (`deeprank` and `lightdock`)

- `deeprank`: Automatically configures `CUDA_VISIBLE_DEVICES` to bind PyTorch graph neural network scoring jobs to allocated GPU hardware.
- `lightdock`: Enables the native `-g` flag for GPU-accelerated glowworm swarm sampling.

---

## 4. CNS CUDA Integration (`cns_solve_CUDA`)

HADDOCK3 can seamlessly integrate with the CUDA-accelerated CNS solver (`cns_solve_CUDA`), accelerating non-bonded energy evaluations and simulated annealing protocols (`rigidbody`, `flexref`, `emref`).

### Building `cns_solve_CUDA`

A build script is provided in the repository to compile the CUDA kernels with automatic GPU architecture detection:

```bash
bash varia/build_cns_cuda.sh
```

The script detects active NVIDIA compute capabilities (`sm_75` for T4, `sm_80` for A100, `sm_86` for A10G, `sm_90` for H100) and compiles the binary to `src/haddock/bin/cns_solve_CUDA`.

### Activating in HADDOCK3

Expose the binary via environment variable or specify it directly:

```bash
export CNS_CUDA_EXEC="/path/to/cns_solve_CUDA"
```

Or set in your configuration file:

```yaml
rigidbody:
  cns_exec: "/path/to/cns_solve_CUDA"
  use_gpu: true
```

### NVIDIA Multi-Process Service (MPS)

When multiplexing multiple concurrent CNS processes onto a single GPU, enable the NVIDIA MPS control daemon to eliminate CUDA context switching latency:

```bash
nvidia-cuda-mps-control -d
```

HADDOCK3's parallel scheduler automatically monitors and interfaces with MPS when present.

---

## 5. High-Performance Computing (HPC) & SLURM Integration

When running in SLURM environments, HADDOCK3 automatically injects GPU resource requests into batch job headers:

```yaml
mode: "slurm"
use_gpu: true
gpus: 2                     # Injects #SBATCH --gres=gpu:2
```

---

## 6. Cloud Testing and Benchmarking with Modal

To run automated test suites and benchmarks on cloud NVIDIA GPUs without local hardware:

```bash
# Install Modal client
pip install 'haddock3[modal]'
modal setup

# Execute unit test suite on an NVIDIA T4 GPU
modal run tests/modal_gpu/modal_runner.py --gpu-type T4 --test-target tests/test_libgpu.py

# Run RMSD matrix benchmark on an NVIDIA A100 GPU
modal run tests/modal_gpu/modal_runner.py --gpu-type A100 --benchmark rmsdmatrix --n-models 2000

# Execute full BM5 benchmark breakdown
modal run tests/modal_gpu/benchmark_bm5.py --target-category rigid --gpu-type A100
```
