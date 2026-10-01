"""Modal cloud GPU backend for HADDOCK3 execution.

This file defines the remote Modal application and container image that
executes GPU-accelerated HADDOCK3 docking jobs on on-demand NVIDIA GPUs.
"""

import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

try:
    import modal
except ImportError:
    modal = None  # type: ignore

# Base repository root for image packaging
_parents = Path(__file__).resolve().parents
REPO_ROOT = _parents[1] if len(_parents) > 1 else Path("/root/haddock3")

if modal is not None:
    app = modal.App("haddock3-gpu-suite")

    # Define the remote GPU container image
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
            "numpy>=1.24.0",
            "scipy>=1.10.0",
            "biopython>=1.80",
            "pandas>=2.0.0",
        )
        .add_local_dir(
            local_path=str(REPO_ROOT),
            remote_path="/root/haddock3",
            ignore=[
                "*.git*",
                "*__pycache__*",
                "*.pytest_cache*",
                "*personal_docs*",
                "*.venv*",
                "*benchmark*",
                "*server/workspace*",
            ],
            copy=True,
        )
        .run_commands(
            "cd /root/haddock3 && pip install --no-build-isolation -e '.[gpu]'"
        )
    )

    @app.function(
        image=haddock_gpu_image,
        gpu="A100",
        timeout=3600,
    )
    def execute_docking_job(
        job_id: str,
        pdb_files: dict[str, str],
        tbl_files: dict[str, str],
        sampling: int = 100,
        refinement: int = 20,
        gpu_device: int = 0,
        gpu_platform: str = "auto",
        mol_params: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        """Execute a full GPU-accelerated docking pipeline on an NVIDIA GPU.

        Args:
            job_id: Unique job identifier.
            pdb_files: Dictionary mapping filename to PDB file text.
            tbl_files: Dictionary mapping filename to restraints text.
            sampling: Number of models to generate at rigid body.
            refinement: Number of top models to refine with flexible annealing.
            gpu_device: GPU device index (e.g. 0).
            gpu_platform: Acceleration backend ('auto', 'cuda').
            mol_params: Optional molecule-specific configuration (e.g. cyclic peptide).

        Returns:
            Dictionary containing structured CAPRI scores, top PDBs, and timing.
        """
        import csv
        import torch

        device_name = (
            torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU"
        )
        work_dir = Path(f"/tmp/haddock_job_{job_id}")
        if work_dir.exists():
            shutil.rmtree(work_dir)
        data_dir = work_dir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        # 1. Write uploaded molecular files
        saved_molecules = []
        for fname, content in pdb_files.items():
            fpath = data_dir / fname
            fpath.write_text(content)
            saved_molecules.append(f"data/{fname}")

        ambig_fname = ""
        for fname, content in tbl_files.items():
            fpath = data_dir / fname
            fpath.write_text(content)
            ambig_fname = f"data/{fname}"

        # 2. Generate HADDOCK3 workflow configuration.
        # NOTE: use_gpu / gpu_devices / gpu_platform in the cfg refer to a CNS
        # CUDA binary that is NOT present in the Modal container. GPU acceleration
        # for FCC/RMSD clustering happens automatically via PyTorch when a GPU is
        # detected. Setting use_gpu=true here breaks CNS and causes 100% job failure.
        cfg_lines = [
            f'run_dir = "run_output"',
            'mode = "local"',
            "ncores = 8",
            f"molecules = {saved_molecules}",
            "",
            "[topoaa]",
            "autohis = false",
        ]
        if mol_params:
            for mol_idx, params in mol_params.items():
                if params.get("cyclic_peptide"):
                    cfg_lines.append("")
                    cfg_lines.append(f"[topoaa.{mol_idx}]")
                    cfg_lines.append("cyclicpept = true")

        # rigidbody stage
        cfg_lines.extend(
            [
                "",
                "[rigidbody]",
                "tolerance = 50",
                f"sampling = {max(sampling, 10)}",
            ]
        )
        if ambig_fname:
            cfg_lines.append(f'ambig_fname = "{ambig_fname}"')
        else:
            # Ab-initio mode: Center-of-mass restraints with screened electrostatics
            # to prevent steric/Coulomb explosion during simulated annealing.
            cfg_lines.append("cmrest = true")
            cfg_lines.append("epsilon = 78")
            cfg_lines.append('dielec = "cdie"')
            cfg_lines.append("randremoval = false")
            cfg_lines.append("w_desolv = 0")

        # flexref stage
        cfg_lines.extend(
            [
                "",
                "[caprieval]",
                "",
                "[seletop]",
                f"select = {max(refinement, 5)}",
                "",
                "[flexref]",
                "tolerance = 50",
            ]
        )
        if ambig_fname:
            cfg_lines.append(f'ambig_fname = "{ambig_fname}"')
        else:
            cfg_lines.append("cmrest = true")
            cfg_lines.append("epsilon = 78")
            cfg_lines.append('dielec = "cdie"')
            cfg_lines.append("randremoval = false")
            cfg_lines.append("w_desolv = 0")

        # emref stage: CNS energy minimisation to clean up remaining clashes
        cfg_lines.extend(
            [
                "",
                "[emref]",
                "tolerance = 50",
            ]
        )
        if ambig_fname:
            cfg_lines.append(f'ambig_fname = "{ambig_fname}"')
        else:
            cfg_lines.append("contactairs = true")
            cfg_lines.append("epsilon = 78")
            cfg_lines.append('dielec = "cdie"')
            cfg_lines.append("randremoval = false")
            cfg_lines.append("w_desolv = 0")

        cfg_lines.extend(
            [
                "",
                "[clustfcc]",
                "min_population = 1",
                "",
                "[rmsdmatrix]",
                "",
                "[seletopclusts]",
                "top_models = 4",
                "",
                "[caprieval]",
                "",
            ]
        )

        cfg_path = work_dir / "workflow.cfg"
        cfg_path.write_text("\n".join(cfg_lines))

        # 3. Execute HADDOCK3
        t0 = time.time()
        cmd = ["haddock3", "workflow.cfg"]
        proc = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            check=False,
        )
        elapsed = round(time.time() - t0, 2)

        out_run_dir = work_dir / "run_output"
        if proc.returncode != 0 or not out_run_dir.exists():
            combined_log = (
                (proc.stdout or "") + "\n" + (proc.stderr or "")
            ).strip()

            diagnostic_details = []
            if out_run_dir.exists():
                # Check for run_output/log
                log_file = out_run_dir / "log"
                if log_file.exists():
                    try:
                        log_txt = log_file.read_text(errors="replace")
                        diagnostic_details.append(f"=== HADDOCK LOG (tail) ===\n{log_txt[-3000:]}")
                    except Exception:
                        pass

                # Check for any .fail files
                for fail_f in sorted(out_run_dir.rglob("*.fail")):
                    try:
                        diagnostic_details.append(f"=== {fail_f.name} ===\n{fail_f.read_text(errors='replace')[:2000]}")
                        break
                    except Exception:
                        pass

                # Check for .cnserr / .cnserr.gz
                import gzip
                for err_f in sorted(out_run_dir.rglob("*.cnserr*")):
                    try:
                        if err_f.name.endswith(".gz"):
                            err_txt = gzip.decompress(err_f.read_bytes()).decode("utf-8", errors="replace")
                        else:
                            err_txt = err_f.read_text(errors="replace")
                        diagnostic_details.append(f"=== {err_f.name} ===\n{err_txt[-3000:]}")
                        break
                    except Exception:
                        pass

                # Check for .out / .out.gz with errors
                for out_f in sorted(out_run_dir.rglob("*.out*"), key=lambda p: p.stat().st_mtime, reverse=True):
                    try:
                        if out_f.name.endswith(".gz"):
                            out_txt = gzip.decompress(out_f.read_bytes()).decode("utf-8", errors="replace")
                        else:
                            out_txt = out_f.read_text(errors="replace")
                        if any(k in out_txt for k in ("%CNS", "ERROR", "^^^^", "ABORT", "error")):
                            diagnostic_details.append(f"=== {out_f.name} ===\n{out_txt[-3000:]}")
                            break
                    except Exception:
                        pass

            error_full = combined_log[-3000:]
            if diagnostic_details:
                error_full += "\n\n" + "\n\n".join(diagnostic_details)

            return {
                "job_id": job_id,
                "status": "failed",
                "gpu_hardware": device_name,
                "execution_seconds": elapsed,
                "returncode": proc.returncode,
                "error_message": error_full[-8000:],
                "clusters": [],
                "best_models": [],
                "pdb_artifacts": {},
            }

        # 4. Parse CAPRI cluster metrics from final step
        clusters = []
        capri_clt_files = sorted(out_run_dir.glob("*caprieval/capri_clt.tsv"))
        if capri_clt_files:
            capri_clt_file = capri_clt_files[-1]
            with open(capri_clt_file, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if not row or row[0].startswith("#") or row[0] == "cluster_rank":
                        continue
                    try:
                        clusters.append(
                            {
                                "cluster_rank": int(row[0]),
                                "cluster_id": int(row[1]),
                                "model_count": int(row[2]),
                                "haddock_score": float(row[4]),
                                "irmsd": float(row[6]),
                                "fnat": float(row[8]),
                                "lrmsd": float(row[10]),
                                "dockq": float(row[12]),
                            }
                        )
                    except (ValueError, IndexError):
                        continue

        # 5. Extract top docked PDB structures (supporting .pdb and .pdb.gz)
        pdb_artifacts = {}
        import gzip

        def read_pdb(p: Path) -> str:
            if p.name.endswith(".gz"):
                return gzip.decompress(p.read_bytes()).decode("utf-8", errors="replace")
            return p.read_text(errors="replace")

        # Priority 1: Top cluster models from seletopclusts
        cluster_pdbs = sorted(out_run_dir.rglob("cluster_*.pdb*"))
        if cluster_pdbs:
            for pdb_path in cluster_pdbs[:10]:
                clean_name = pdb_path.name.removesuffix(".gz")
                pdb_artifacts[clean_name] = read_pdb(pdb_path)
        else:
            # Priority 2: Any refined models (emref or flexref)
            refined_pdbs = sorted(out_run_dir.rglob("*emref*/*.pdb*")) or sorted(out_run_dir.rglob("*flexref*/*.pdb*"))
            for pdb_path in refined_pdbs[:10]:
                if not any(pdb_path.name.endswith(ext) for ext in (".cnserr", ".cnserr.gz", ".out", ".out.gz")):
                    clean_name = pdb_path.name.removesuffix(".gz")
                    pdb_artifacts[clean_name] = read_pdb(pdb_path)

        # Fallback: Any PDB files in run_output
        if not pdb_artifacts:
            for pdb_path in sorted(out_run_dir.rglob("*.pdb*"))[:10]:
                if not any(pdb_path.name.endswith(ext) for ext in (".cnserr", ".cnserr.gz", ".out", ".out.gz")):
                    clean_name = pdb_path.name.removesuffix(".gz")
                    pdb_artifacts[clean_name] = read_pdb(pdb_path)

        return {
            "job_id": job_id,
            "status": "completed",
            "gpu_hardware": device_name,
            "execution_seconds": elapsed,
            "returncode": 0,
            "clusters": clusters,
            "best_models": list(pdb_artifacts.keys()),
            "pdb_artifacts": pdb_artifacts,
            "log_tail": proc.stdout[-2000:] if proc.stdout else "",
        }
