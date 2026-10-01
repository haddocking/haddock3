"""System health and GPU infrastructure inspection endpoints."""

from datetime import datetime, timezone
from fastapi import APIRouter

from server.models import GpuInfoResponse

router = APIRouter(prefix="/system", tags=["System"])


@router.get(
    "/health",
    summary="Service health check",
    description="Returns status of the FastAPI control plane running on Google Cloud Run.",
)
async def health_check() -> dict[str, str]:
    """Provide health probe response for Cloud Run liveness/readiness checks."""
    return {
        "status": "healthy",
        "service": "haddock3-gpu-control-plane",
        "version": "3.0.0-gpu",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


@router.get(
    "/gpu",
    response_model=GpuInfoResponse,
    summary="Get available GPU compute infrastructure",
    description="Inspect supported cloud GPU accelerators (A100, H100, A10G, T4) available on Modal.",
)
async def get_gpu_info() -> GpuInfoResponse:
    """Return backend capabilities, default GPU architecture, and deployment environment."""
    return GpuInfoResponse(
        backend="Modal Cloud GPU Serverless Engine",
        supported_gpus=["T4", "A10G", "A100", "H100"],
        default_gpu="A100",
        service_status="healthy",
        cloud_run_environment=True,
    )
