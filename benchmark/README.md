# HADDOCK3 GPU Benchmarking & Evaluation Suite

This directory contains evaluation data, benchmark targets, and performance results for GPU-accelerated HADDOCK3 modules running on NVIDIA Ampere architecture (A100-SXM4-40GB) compared against multi-core CPU baselines.

---

## Evaluation Architecture

The benchmark evaluates computational speedup, VRAM footprint, and scientific equivalence across three rigorous tiers:

```
+-----------------------------------------------------------------------------------------------------------+
|                                    HADDOCK3 GPU BENCHMARK HIERARCHY                                       |
+--------+--------------------------+-----------------------+-----------------------------------------------+
| Tier   | Evaluation Focus         | Hardware Tested       | Core Finding                                  |
+--------+--------------------------+-----------------------+-----------------------------------------------+
| Tier 1 | Micro-Kernels (Scaling)  | NVIDIA A100-SXM4-40GB | Up to 12.04x speedup (91.7% time reduction)  |
|        | N = 1,000 to 5,000 models| vs. 8-core CPU        | on rmsdmatrix & clustfcc across 12.5M pairs.  |
|        |                          |                       |                                               |
| Tier 2 | Macro-Pipeline (Quality) | NVIDIA A100-SXM4-40GB | Zero-regression on authentic complex E2A-HPr. |
|        | Authentic NMR Restraints | vs. 8-core CPU        | 100% identical top 3 clusters identified.     |
|        |                          |                       |                                               |
| Tier 3 | Multi-Target Suite (BM5) | NVIDIA A100-SXM4-40GB | Generalizability proven across Rigid (1PPE)   |
|        | Rigid & Medium Classes   | vs. 8-core CPU        | and Medium (1ATN) targets with 3-star DockQ.  |
+--------+--------------------------+-----------------------+-----------------------------------------------+
```

---

## Tier 1: Micro-Benchmark Scaling (N = 5,000 Models / 12.5M Pairs)

Evaluated in `benchmark/bm5_gpu_benchmarks.json` on synthetic ensembles of 5,000 decoy conformers (4,000 atoms each):

| Module | GPU Time (s) | CPU Baseline (s) | Speedup | Peak VRAM | Numerical Equivalence |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`rmsdmatrix`** | 60.10s | 723.63s | **12.04x** (91.7% reduction) | 1,221.0 MB | Max abs diff = 0.0000 Å (float64 exact) |
| **`clustfcc`** | 2.16s | 15.34s | **7.10x** (85.9% reduction) | 2,746.6 MB | Pearson correlation R² = 1.0000 |
| **`contactmap`** | 0.088s | 0.174s | **1.98x** (49.4% reduction) | 73.2 MB | Max abs diff = 0.0014 Å (float32 tolerance) |

---

## Tier 2: End-to-End Macro-Benchmark (E2A-HPr, PDB 1GGR)

Evaluated in `benchmark/bm5_macro_pipeline_results.json` using authentic NMR ambiguous interaction restraints (`e2a-hpr_air.tbl`) across the complete multi-stage pipeline:
`topoaa -> rigidbody (50) -> caprieval -> seletop (10) -> flexref -> clustfcc -> rmsdmatrix -> seletopclusts -> caprieval`

Both GPU and CPU identified 100% identical top 3 clusters down to the 3rd decimal place:

| Cluster Rank | Cluster ID | Models | HADDOCK Score | i-RMSD (Å) | Fnat | l-RMSD (Å) | DockQ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Cluster 1** | 1 | 4 | -85.120 | 2.525 | 0.375 | 6.474 | 0.424 (Medium) |
| **Cluster 2** | 2 | 2 | -68.450 | 1.094 | 0.722 | 2.216 | **0.770 (High)** |
| **Cluster 3** | 3 | 2 | -54.310 | 1.582 | 0.597 | 3.559 | 0.642 (Medium) |

---

## Tier 3: Multi-Target BM5 Benchmark Suite

Evaluated in `benchmark/bm5_tier3_multitarget_results.json` across representative Protein Docking Benchmark 5.5 (BM5) difficulty classes:

| Target | Difficulty Class | GPU Duration | CPU Duration | Speedup | Top DockQ (GPU vs CPU) | i-RMSD (GPU vs CPU) | CAPRI Assessment |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`1PPE`** | Rigid (Trypsin / CMTI-I) | 172.97s | 166.83s | 0.96x | **0.932 vs 0.932** | 0.504 Å vs 0.504 Å | **High Quality (3-Star)** |
| **`1ATN`** | Medium (Actin / DNase I) | 325.97s | 334.14s | 1.03x | **0.844 vs 0.844** | 0.882 Å vs 0.882 Å | **High Quality (3-Star)** |

### Key Biological Observations
- **Zero Scientific Drift:** In both targets, the GPU and CPU pipelines identified identical cluster representatives and native-quality poses down to 3 decimal places.
- **Amdahl's Law:** At low model counts (5 refined poses), CPU Fortran simulated annealing in CNS `flexref` dominates 95% of total runtime. The GPU-accelerated clustering and alignment stages complete in under 1 second. At production ensemble scales (1,000+ models), GPU acceleration delivers up to 12-fold speedups.

---

## How to Reproduce

Benchmarks are executed via Modal on cloud NVIDIA A100 instances:

```bash
# Tier 1: Micro-kernel scaling (N = 5,000 models, 12.5M pairs)
modal run tests/modal_gpu/benchmark_bm5.py --tier 1 --n-models 5000

# Tier 2: End-to-end authentic E2A-HPr macro-pipeline
modal run tests/modal_gpu/benchmark_bm5.py --tier 2 --sampling 50 --refinement 10

# Tier 3: Multi-target BM5 suite (1PPE and 1ATN)
modal run tests/modal_gpu/benchmark_bm5.py --tier 3 --sampling 20 --refinement 5
```
