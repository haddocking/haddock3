"""Modal cloud GPU execution and benchmarking harness for HADDOCK3.

Allows running HADDOCK3 GPU tests, refinement, and matrix benchmarks on
real NVIDIA GPUs (T4, A10G, A100, H100) via Modal.

Usage:
    modal run tests/modal_gpu/modal_runner.py --gpu-type T4 --test-target tests/test_libgpu.py
    modal run tests/modal_gpu/modal_runner.py --gpu-type A10G --benchmark rmsdmatrix
"""

import subprocess
import sys
import time
from pathlib import Path

try:
    import modal
except ImportError:
    modal = None  # type: ignore

# Base workspace path for mounting local source
REPO_ROOT = Path(__file__).resolve().parents[2]

if modal is not None:
    app = modal.App("haddock3-gpu-suite")

    # Define container image with CUDA, PyTorch, OpenMM, and HADDOCK3
    haddock_gpu_image = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install(
            "git",
            "build-essential",
            "gfortran",
            "tcsh",
            "libopenmpi-dev",
            "openmpi-bin",
            "curl",
        )
        .pip_install(
            "torch>=2.0.0",
            "openmm>=8.0.0",
            "pdbfixer",
            "pytest>=8.0.0",
            "pytest-mock",
            "pytest-cov",
            "numpy>=1.24.0",
            "scipy>=1.10.0",
            "biopython>=1.80",
        )
        .add_local_dir(
            local_path=str(REPO_ROOT),
            remote_path="/root/haddock3",
            ignore=["*.git*", "*__pycache__*", "*.pytest_cache*", "*personal_docs*"],
        )
        .run_commands(
            "cd /root/haddock3 && pip install --no-build-isolation -e '.[gpu]'"
        )
    )

    @app.function(
        image=haddock_gpu_image,
        gpu="T4",
        timeout=1800,
    )
    def run_tests_on_gpu(
        test_path: str = "tests/test_libgpu.py",
    ) -> dict[str, str | int]:
        """Execute pytest test suite on remote NVIDIA GPU.

        Args:
            test_path: Relative path to test file or directory.

        Returns:
            dict containing returncode, stdout, stderr, and GPU device details.
        """
        import torch

        gpu_name = (
            torch.cuda.get_device_name(0) if torch.cuda.is_available() else "None"
        )
        cuda_ver = torch.version.cuda or "None"

        cmd = [
            "pytest",
            "-v",
            f"/root/haddock3/{test_path}",
        ]

        start_time = time.time()
        result = subprocess.run(
            cmd,
            cwd="/root/haddock3",
            capture_output=True,
            text=True,
            check=False,
        )
        elapsed = time.time() - start_time

        return {
            "gpu_name": gpu_name,
            "cuda_version": cuda_ver,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "elapsed_seconds": round(elapsed, 2),
        }

    @app.function(
        image=haddock_gpu_image,
        gpu="A100",
        timeout=3600,
    )
    def benchmark_kernel(
        kernel_name: str = "rmsdmatrix",
        n_models: int = 1000,
    ) -> dict[str, float | str | int]:
        """Benchmark GPU vs CPU performance for compute-intensive modules.

        Args:
            kernel_name: Name of module/kernel to benchmark ('rmsdmatrix', 'clustfcc', 'openmm').
            n_models: Number of conformational models to simulate.

        Returns:
            dict with speedup factor, GPU execution time, CPU execution time.
        """
        import time

        import numpy as np
        import torch

        if not torch.cuda.is_available():
            raise RuntimeError("No GPU device detected in benchmark runner.")

        device_name = torch.cuda.get_device_name(0)

        if kernel_name == "rmsdmatrix":
            from haddock.libs.libalign_gpu import compute_rmsd_matrix

            coords = np.random.randn(n_models, 250, 3).astype(np.float64)
            _ = compute_rmsd_matrix(coords[:10], use_gpu=True, device="cuda")

            t0 = time.perf_counter()
            _, _, gpu_rmsds = compute_rmsd_matrix(
                coords, use_gpu=True, device="cuda"
            )
            t_gpu = time.perf_counter() - t0

            t0 = time.perf_counter()
            _, _, _ = compute_rmsd_matrix(coords, use_gpu=False, device="cpu")
            t_cpu = time.perf_counter() - t0

            speedup = t_cpu / max(t_gpu, 1e-6)
            return {
                "kernel": "rmsdmatrix",
                "device": device_name,
                "n_models": n_models,
                "n_pairs": len(gpu_rmsds),
                "gpu_seconds": round(t_gpu, 4),
                "cpu_seconds": round(t_cpu, 4),
                "speedup": round(speedup, 2),
            }

        elif kernel_name == "clustfcc":
            from haddock.libs.libfcc import calculate_pairwise_matrix
            from haddock.libs.libfcc_gpu import calculate_pairwise_matrix_gpu

            pool = list(range(2000))
            contacts = [
                set(
                    np.random.choice(
                        pool, size=np.random.randint(40, 120), replace=False
                    )
                )
                for _ in range(n_models)
            ]

            _ = calculate_pairwise_matrix_gpu(contacts[:10], device="cuda")

            t0 = time.perf_counter()
            _ = calculate_pairwise_matrix_gpu(contacts, device="cuda")
            t_gpu = time.perf_counter() - t0

            cpu_sample = contacts[: min(n_models, 300)]
            t0 = time.perf_counter()
            _ = list(calculate_pairwise_matrix(cpu_sample, ignore_chain=False))
            t_cpu_sample = time.perf_counter() - t0

            total_pairs = n_models * (n_models - 1) // 2
            sample_pairs = len(cpu_sample) * (len(cpu_sample) - 1) // 2
            t_cpu = t_cpu_sample * (total_pairs / max(sample_pairs, 1))

            speedup = t_cpu / max(t_gpu, 1e-6)
            return {
                "kernel": "clustfcc",
                "device": device_name,
                "n_models": n_models,
                "n_pairs": total_pairs,
                "gpu_seconds": round(t_gpu, 4),
                "cpu_extrapolated_seconds": round(t_cpu, 4),
                "speedup": round(speedup, 2),
            }

        elif kernel_name == "contactmap":
            from scipy.spatial.distance import pdist, squareform

            from haddock.modules.analysis.contactmap.contmap import (
                compute_distance_matrix,
            )

            n_atoms = n_models * 10
            coords = np.random.randn(n_atoms, 3).tolist()

            _ = compute_distance_matrix(
                coords[:100], use_gpu=True, device="cuda"
            )

            t0 = time.perf_counter()
            _ = compute_distance_matrix(coords, use_gpu=True, device="cuda")
            t_gpu = time.perf_counter() - t0

            t0 = time.perf_counter()
            _ = squareform(pdist(coords))
            t_cpu = time.perf_counter() - t0

            speedup = t_cpu / max(t_gpu, 1e-6)
            return {
                "kernel": "contactmap",
                "device": device_name,
                "n_atoms": n_atoms,
                "gpu_seconds": round(t_gpu, 4),
                "cpu_seconds": round(t_cpu, 4),
                "speedup": round(speedup, 2),
            }

        raise ValueError(f"Unknown kernel: {kernel_name}")

    @app.local_entrypoint()
    def main(
        gpu_type: str = "T4",
        test_target: str = "tests/test_libgpu.py",
        benchmark: str = "",
        n_models: int = 500,
    ) -> None:
        """Local CLI entrypoint for modal run."""
        print(f"=== HADDOCK3 Cloud GPU Harness (Target GPU: {gpu_type}) ===")
        if benchmark:
            print(f"Running benchmark '{benchmark}' with {n_models} models on A100...")
            res = benchmark_kernel.remote(kernel_name=benchmark, n_models=n_models)
            print("Benchmark result:", res)
        else:
            print(f"Running test suite '{test_target}' on {gpu_type}...")
            res = run_tests_on_gpu.remote(test_path=test_target)
            print(f"GPU: {res['gpu_name']} (CUDA {res['cuda_version']})")
            print(
                f"Exit code: {res['returncode']} (Duration: {res['elapsed_seconds']}s)"
            )
            print("\n--- Output ---")
            print(res["stdout"])
            if res["stderr"]:
                print("\n--- Errors ---")
                print(res["stderr"])
            if res["returncode"] != 0:
                sys.exit(res["returncode"])
