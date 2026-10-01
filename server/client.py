"""Asynchronous Modal dispatch client and job registry."""

import asyncio
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from server.config import settings
from server.models import CapriClusterMetric, JobResultsResponse, JobStatus, JobStatusResponse

try:
    import modal
    from server.modal_backend import execute_docking_job
except ImportError:
    modal = None
    execute_docking_job = None


class JobRecord:
    """In-memory metadata record for a submitted docking job."""

    def __init__(
        self,
        job_id: str,
        gpu_type: str,
        job_name: Optional[str] = None,
        modal_call_id: Optional[str] = None,
    ):
        self.job_id = job_id
        self.gpu_type = gpu_type
        self.job_name = job_name
        self.modal_call_id = modal_call_id
        self.status: JobStatus = JobStatus.QUEUED
        self.created_at = datetime.now(timezone.utc).isoformat()
        self.started_at: Optional[float] = None
        self.finished_at: Optional[float] = None
        self.current_stage: Optional[str] = "00_topoaa"
        self.completed_stages: list[str] = []
        self.error_message: Optional[str] = None
        self.results: Optional[dict[str, Any]] = None
        self.call_handle: Any = None


class ModalClientManager:
    """Manages async dispatch, polling, and results caching for Modal jobs."""

    def __init__(self):
        self.jobs: dict[str, JobRecord] = {}
        settings.WORKSPACE_DIR.mkdir(parents=True, exist_ok=True)
        self._load_persisted_jobs()

    def _load_persisted_jobs(self) -> None:
        """Scan workspace directory and rehydrate known completed or cached jobs."""
        import json
        for job_dir in settings.WORKSPACE_DIR.iterdir():
            if not job_dir.is_dir():
                continue
            state_file = job_dir / "state.json"
            if state_file.exists():
                try:
                    data = json.loads(state_file.read_text())
                    job_id = data.get("job_id", job_dir.name)
                    record = JobRecord(
                        job_id=job_id,
                        gpu_type=data.get("gpu_type", "H100"),
                        job_name=data.get("job_name"),
                        modal_call_id=data.get("modal_call_id"),
                    )
                    record.status = JobStatus(data.get("status", "completed"))
                    record.current_stage = data.get("current_stage", "completed")
                    record.completed_stages = data.get("completed_stages", [])
                    record.started_at = data.get("started_at", 0)
                    record.finished_at = data.get("started_at", 0) + data.get("execution_seconds", 0)
                    results_file = job_dir / "results.json"
                    if results_file.exists():
                        record.results = json.loads(results_file.read_text())
                    self.jobs[job_id] = record
                except Exception:
                    pass

    async def submit_job(
        self,
        job_id: str,
        pdb_files: dict[str, str],
        tbl_files: dict[str, str],
        sampling: int = 100,
        refinement: int = 20,
        gpu_type: str = "A100",
        gpu_device: int = 0,
        gpu_platform: str = "auto",
        job_name: Optional[str] = None,
        mol_params: Optional[dict[str, Any]] = None,
    ) -> JobRecord:
        """Asynchronously dispatch a docking calculation to Modal."""
        record = JobRecord(job_id=job_id, gpu_type=gpu_type, job_name=job_name)
        record.started_at = time.time()
        record.status = JobStatus.RUNNING

        if modal is not None:
            try:
                # 1. Lookup the deployed remote function from Modal
                try:
                    fn = modal.Function.from_name(settings.MODAL_APP_NAME, "execute_docking_job")
                except Exception:
                    fn = execute_docking_job

                # 2. Dispatch asynchronously using Modal's native async spawn
                if hasattr(fn, "spawn") and hasattr(fn.spawn, "aio"):
                    call = await fn.spawn.aio(
                        job_id=job_id,
                        pdb_files=pdb_files,
                        tbl_files=tbl_files,
                        sampling=sampling,
                        refinement=refinement,
                        gpu_device=gpu_device,
                        gpu_platform=gpu_platform,
                        mol_params=mol_params,
                    )
                elif hasattr(fn, "spawn"):
                    call = fn.spawn(
                        job_id=job_id,
                        pdb_files=pdb_files,
                        tbl_files=tbl_files,
                        sampling=sampling,
                        refinement=refinement,
                        gpu_device=gpu_device,
                        gpu_platform=gpu_platform,
                        mol_params=mol_params,
                    )
                else:
                    raise RuntimeError("No executable Modal function available.")

                record.modal_call_id = getattr(call, "object_id", str(call))
                record.call_handle = call
            except Exception as e:
                # If running locally or without credentials, register error
                record.status = JobStatus.FAILED
                record.error_message = f"Modal dispatch error: {str(e)}"
        else:
            record.status = JobStatus.FAILED
            record.error_message = "Modal SDK is not installed or configured."

        self.jobs[job_id] = record
        return record

    async def get_status(self, job_id: str) -> Optional[JobStatusResponse]:
        """Poll the asynchronous state of a running job."""
        record = self.jobs.get(job_id)
        if not record:
            return None

        # Check Modal call status if still running
        if record.status == JobStatus.RUNNING and record.call_handle:
            try:
                # Non-blocking get
                if hasattr(record.call_handle, "get") and hasattr(record.call_handle.get, "aio"):
                    res = await record.call_handle.get.aio(timeout=0)
                elif hasattr(record.call_handle, "get"):
                    res = record.call_handle.get(timeout=0)
                else:
                    res = None

                if res is not None:
                    record.results = res
                    if res.get("status") == "completed":
                        record.status = JobStatus.COMPLETED
                        record.finished_at = time.time()
                        record.current_stage = "completed"
                        record.completed_stages = [
                            "00_topoaa",
                            "01_rigidbody",
                            "02_caprieval",
                            "03_seletop",
                            "04_flexref",
                            "05_emref",
                            "06_clustfcc",
                            "07_rmsdmatrix",
                            "08_seletopclusts",
                            "09_caprieval",
                        ]
                        self._persist_artifacts(job_id, res)
                    else:
                        record.status = JobStatus.FAILED
                        record.error_message = res.get("error_message", "Unknown error")
            except TimeoutError:
                # Still running
                pass
            except Exception as e:
                record.status = JobStatus.FAILED
                record.error_message = str(e)

        if record.status == JobStatus.COMPLETED and not record.completed_stages:
            record.current_stage = "completed"
            record.completed_stages = [
                "00_topoaa",
                "01_rigidbody",
                "02_caprieval",
                "03_seletop",
                "04_flexref",
                "05_emref",
                "06_clustfcc",
                "07_rmsdmatrix",
                "08_seletopclusts",
                "09_caprieval",
            ]

        if record.results and "execution_seconds" in record.results:
            elapsed = round(float(record.results["execution_seconds"]), 2)
        else:
            elapsed = (
                round((record.finished_at or time.time()) - (record.started_at or time.time()), 2)
            )

        return JobStatusResponse(
            job_id=record.job_id,
            status=record.status,
            gpu_type=record.gpu_type,
            current_stage=record.current_stage,
            completed_stages=record.completed_stages,
            elapsed_seconds=elapsed,
            error_message=record.error_message,
        )

    async def get_results(self, job_id: str) -> Optional[JobResultsResponse]:
        """Fetch structured results for a completed job."""
        record = self.jobs.get(job_id)
        if not record or record.status != JobStatus.COMPLETED:
            return None

        res = record.results or {}
        clusters = [
            CapriClusterMetric(
                cluster_rank=c["cluster_rank"],
                cluster_id=c["cluster_id"],
                model_count=c["model_count"],
                haddock_score=c["haddock_score"],
                fnat=c["fnat"],
                irmsd=c["irmsd"],
                lrmsd=c["lrmsd"],
                dockq=c["dockq"],
            )
            for c in res.get("clusters", [])
        ]

        best_models = res.get("best_models", [])
        download_urls = {
            model_name: f"/api/v1/results/{job_id}/download/{model_name}"
            for model_name in best_models
        }

        if "execution_seconds" in res:
            elapsed = round(float(res["execution_seconds"]), 2)
        else:
            elapsed = (
                round((record.finished_at or time.time()) - (record.started_at or time.time()), 2)
            )

        return JobResultsResponse(
            job_id=job_id,
            status=record.status,
            gpu_type=record.gpu_type,
            execution_seconds=elapsed,
            clusters=clusters,
            best_models=best_models,
            download_urls=download_urls,
        )

    async def cancel_job(self, job_id: str) -> bool:
        """Cancel a running Modal job."""
        record = self.jobs.get(job_id)
        if not record:
            return False

        if record.call_handle:
            try:
                record.call_handle.cancel()
            except Exception:
                pass

        record.status = JobStatus.CANCELLED
        return True

    def _persist_artifacts(self, job_id: str, res: dict[str, Any]) -> None:
        """Persist downloaded PDBs to local job workspace for streaming downloads."""
        import json
        job_dir = settings.WORKSPACE_DIR / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        artifacts = res.get("pdb_artifacts", {})
        for fname, content in artifacts.items():
            fpath = job_dir / fname
            fpath.write_text(content)

        (job_dir / "results.json").write_text(json.dumps(res, indent=2))

        record = self.jobs.get(job_id)
        if record:
            state = {
                "job_id": job_id,
                "gpu_type": record.gpu_type,
                "job_name": record.job_name,
                "modal_call_id": record.modal_call_id,
                "status": record.status.value,
                "current_stage": record.current_stage,
                "completed_stages": record.completed_stages,
                "started_at": record.started_at,
                "execution_seconds": res.get("execution_seconds"),
            }
            (job_dir / "state.json").write_text(json.dumps(state, indent=2))


# Global singleton manager instance
modal_manager = ModalClientManager()
