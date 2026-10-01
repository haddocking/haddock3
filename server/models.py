"""Data models and schemas for HADDOCK3 FastAPI Cloud Service."""

from enum import Enum
from typing import Any, Optional
from pydantic import BaseModel, Field


class JobStatus(str, Enum):
    """Execution status of a docking job."""

    QUEUED = "queued"
    INITIALIZING = "initializing"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class GpuType(str, Enum):
    """Supported cloud GPU accelerators on Modal."""

    T4 = "T4"
    A10G = "A10G"
    A100 = "A100"
    H100 = "H100"


class JobSubmissionResponse(BaseModel):
    """Response returned immediately upon successful job submission."""

    job_id: str = Field(..., description="Unique UUID identifying the docking job.")
    status: JobStatus = Field(JobStatus.QUEUED, description="Current job lifecycle state.")
    gpu_type: str = Field(..., description="Target NVIDIA GPU allocated on Modal.")
    modal_call_id: Optional[str] = Field(None, description="Remote Modal FunctionCall ID.")
    status_url: str = Field(..., description="Endpoint to poll for execution progress.")
    results_url: str = Field(..., description="Endpoint to fetch final scientific results.")
    created_at: str = Field(..., description="ISO 8601 UTC submission timestamp.")
    error_message: Optional[str] = Field(None, description="Detailed error description if submission failed.")


class JobStatusResponse(BaseModel):
    """Detailed real-time execution progress of a docking job."""

    job_id: str
    status: JobStatus
    gpu_type: str
    current_stage: Optional[str] = Field(
        None, description="Active HADDOCK3 workflow module (e.g. 05_clustfcc)."
    )
    completed_stages: list[str] = Field(
        default_factory=list, description="List of finished pipeline stages."
    )
    total_stages: int = Field(10, description="Total stages in the docking workflow.")
    elapsed_seconds: float = Field(0.0, description="Elapsed wall-clock execution time.")
    error_message: Optional[str] = Field(None, description="Error detail if job failed.")


class CapriClusterMetric(BaseModel):
    """CAPRI scientific evaluation metrics for a docked cluster."""

    cluster_rank: int = Field(..., description="Ranking ordered by HADDOCK score.")
    cluster_id: int = Field(..., description="Cluster identifier from clustfcc.")
    model_count: int = Field(..., description="Number of models assigned to this cluster.")
    haddock_score: float = Field(..., description="Weighted HADDOCK docking score (lower is better).")
    fnat: float = Field(..., description="Fraction of native contacts preserved (0.0 to 1.0).")
    irmsd: float = Field(..., description="Interface backbone RMSD to native in Angstroms.")
    lrmsd: float = Field(..., description="Ligand backbone RMSD to native in Angstroms.")
    dockq: float = Field(..., description="Standard DockQ overall quality score (0.0 to 1.0).")


class JobResultsResponse(BaseModel):
    """Final docked structures, CAPRI rankings, and downloadable artifacts."""

    job_id: str
    status: JobStatus
    gpu_type: str
    execution_seconds: float
    clusters: list[CapriClusterMetric] = Field(
        default_factory=list, description="CAPRI cluster metrics."
    )
    best_models: list[str] = Field(
        default_factory=list, description="Filenames of top docked PDB structures."
    )
    download_urls: dict[str, str] = Field(
        default_factory=dict, description="Named links for structure downloads and reports."
    )


class GpuInfoResponse(BaseModel):
    """System information regarding cloud GPU backends."""

    backend: str = "Modal Cloud GPU Serverless Engine"
    supported_gpus: list[str] = ["T4", "A10G", "A100", "H100"]
    default_gpu: str = "A100"
    service_status: str = "healthy"
    cloud_run_environment: bool = True
