"""Results retrieval and structural artifact download endpoints."""

from pathlib import Path
from fastapi import APIRouter, HTTPException, status
from fastapi.responses import FileResponse

from server.client import modal_manager
from server.config import settings
from server.models import JobResultsResponse, JobStatus

router = APIRouter(prefix="/results", tags=["Results"])


@router.get(
    "/{job_id}",
    response_model=JobResultsResponse,
    summary="Retrieve CAPRI rankings and top docked structures for a completed job",
)
async def get_docking_results(job_id: str) -> JobResultsResponse:
    """Fetch scientific metrics (DockQ, i-RMSD, l-RMSD, fnat, HADDOCK-score) and model filenames."""
    record = modal_manager.jobs.get(job_id)
    if not record:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job '{job_id}' not found in registry.",
        )

    if record.status in (JobStatus.QUEUED, JobStatus.INITIALIZING, JobStatus.RUNNING):
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=(
                f"Job '{job_id}' is currently {record.status.value}. "
                f"Please poll /api/v1/docking/{job_id}/status until completion."
            ),
        )

    if record.status == JobStatus.FAILED:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Job '{job_id}' failed: {record.error_message or 'Unknown error'}",
        )

    if record.status == JobStatus.CANCELLED:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Job '{job_id}' was cancelled by user.",
        )

    results = await modal_manager.get_results(job_id)
    if not results:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Results for job '{job_id}' are not yet available.",
        )

    return results


@router.get(
    "/{job_id}/download/{filename}",
    response_class=FileResponse,
    summary="Download docked PDB structure or scientific report artifact",
)
async def download_structure_file(job_id: str, filename: str) -> FileResponse:
    """Stream a specific docked PDB structure or report file for a given job."""
    # Sanitize filename against directory traversal attacks
    safe_filename = Path(filename).name
    if safe_filename != filename or ".." in filename:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid filename specified.",
        )

    file_path = settings.WORKSPACE_DIR / job_id / safe_filename
    if not file_path.is_file():
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"File '{safe_filename}' not found for job '{job_id}'.",
        )

    media_type = "chemical/x-pdb" if safe_filename.endswith(".pdb") else "text/plain"
    return FileResponse(
        path=file_path,
        media_type=media_type,
        filename=safe_filename,
    )
