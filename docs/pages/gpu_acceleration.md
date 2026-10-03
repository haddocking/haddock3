# GPU Acceleration in HADDOCK3

HADDOCK3 includes optional, zero-regression GPU acceleration designed to accelerate pairwise $O(N^2)$ post-docking analysis and clustering (RMSD matrices, Fraction of Common Contacts, contact maps) and molecular dynamics refinement.

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
- Achieves up to ~12x speedup over CPU execution for large ensembles ($N \ge 1,000$ structures) with double-precision numerical equivalence (maximum absolute difference = 0.0000 Å vs CPU baseline).

### 3.2 Fraction of Common Contacts Clustering (`clustfcc`)

The Fraction of Common Contacts (FCC) algorithm calculates the residue-residue contact overlap between all pairs of decoy complexes. In pure Python, evaluating $N \times (N - 1) / 2$ set intersections is a major bottleneck for large ensembles.

With `use_gpu: true`:
- Contacts are encoded as a binary occurrence matrix $A$ of shape $(N, \text{total\_unique\_contacts})$.
- Pairwise intersection counts are evaluated simultaneously using matrix multiplication: $M = A \times A^T$.
- Accelerated using GPU Tensor Cores via PyTorch, achieving up to ~7x speedup over the serial CPU implementation with exact cluster preservation ($R^2 = 1.0000$).

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

## 4. Scope and Workflow Performance (Amdahl's Law)

In standard HADDOCK3 docking workflows, sampling and simulated annealing stages (`rigidbody`, `flexref`, `mdref`) rely on the Crystallography & NMR System (CNS) Fortran executable, which runs on CPU. Because CNS dominates 85% to >95% of total wall-clock runtime in a standard docking pipeline, overall end-to-end workflow runtimes remain primarily CPU-bound in accordance with Amdahl's Law.

The GPU acceleration introduced here specifically addresses the memory footprint and execution bottleneck during **post-docking ensemble analysis and clustering**:
- For standard small runs ($N \le 200$), CPU clustering completes in seconds, so GPU overhead yields near parity (~1.0x).
- For production ensembles ($N = 1,000$ to $10,000+$ models), pairwise $O(N^2)$ matrix evaluations otherwise become severe CPU/RAM bottlenecks. GPU execution reduces these analysis stages from tens of minutes to seconds.

---

## 5. High-Performance Computing (HPC) & SLURM Integration

When running in SLURM environments, HADDOCK3 automatically injects GPU resource requests into batch job headers:

```yaml
mode: "slurm"
use_gpu: true
gpus: 2                     # Injects #SBATCH --gres=gpu:2
```
