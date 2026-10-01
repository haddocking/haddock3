"""Docking submission, monitoring, and lifecycle endpoints."""

import uuid
from typing import Optional
from fastapi import APIRouter, File, Form, HTTPException, UploadFile, status

from server.client import modal_manager
from server.config import settings
from server.models import (
    GpuType,
    JobStatus,
    JobStatusResponse,
    JobSubmissionResponse,
)

router = APIRouter(prefix="/docking", tags=["Docking"])


@router.post(
    "/submit",
    response_model=JobSubmissionResponse,
    status_code=status.HTTP_202_ACCEPTED,
    summary="Submit a molecular docking calculation to Modal GPU",
    description=(
        "Submits 2 or more PDB structures and optional restraint files (.tbl/.act) "
        "for accelerated docking on dedicated NVIDIA cloud GPUs via Modal."
    ),
)
async def submit_docking_job(
    molecules: list[UploadFile] = File(
        ...,
        description="At least 2 PDB structures (e.g. receptor and ligand) to dock.",
    ),
    restraints: Optional[UploadFile] = File(
        None,
        description="Optional Ambiguous Interaction Restraint (.tbl or .act) file.",
    ),
    sampling: int = Form(
        100,
        ge=1,
        le=1000,
        description="Number of rigid-body models to generate in sampling stage.",
    ),
    refinement: int = Form(
        20,
        ge=1,
        le=500,
        description="Number of models to refine flexibly with OpenMM GPU.",
    ),
    gpu_type: GpuType = Form(
        GpuType.A100,
        description="Target cloud GPU accelerator on Modal (default: A100).",
    ),
    gpu_device: int = Form(
        0,
        ge=0,
        description="CUDA device index for compute execution.",
    ),
    gpu_platform: str = Form(
        "auto",
        description="OpenMM acceleration platform ('auto', 'CUDA', or 'OpenCL').",
    ),
    job_name: Optional[str] = Form(
        None,
        max_length=120,
        description="Optional descriptive title for the docking run.",
    ),
) -> JobSubmissionResponse:
    """Validate incoming molecular files and dispatch execution to Modal GPUs asynchronously."""
    if len(molecules) < 2:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="At least 2 molecular structures (e.g. receptor and ligand) are required for docking.",
        )

    max_bytes = settings.MAX_UPLOAD_SIZE_MB * 1024 * 1024
    pdb_files: dict[str, str] = {}

    for mol in molecules:
        fname = mol.filename or f"molecule_{len(pdb_files) + 1}.pdb"
        if not (fname.endswith(".pdb") or fname.endswith(".ent") or fname.endswith(".cif")):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Unsupported file format for '{fname}'. Must be .pdb, .ent, or .cif.",
            )

        content_bytes = await mol.read()
        if len(content_bytes) > max_bytes:
            raise HTTPException(
                status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
                detail=f"File '{fname}' exceeds the {settings.MAX_UPLOAD_SIZE_MB}MB size limit.",
            )
        if len(content_bytes) == 0:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"File '{fname}' is empty.",
            )

        try:
            pdb_files[fname] = content_bytes.decode("utf-8")
        except UnicodeDecodeError:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"File '{fname}' must be valid UTF-8 text.",
            )

    tbl_files: dict[str, str] = {}
    if restraints is not None:
        tbl_name = restraints.filename or "ambig.tbl"
        tbl_bytes = await restraints.read()
        if len(tbl_bytes) > max_bytes:
            raise HTTPException(
                status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
                detail=f"Restraints file '{tbl_name}' exceeds the {settings.MAX_UPLOAD_SIZE_MB}MB limit.",
            )
        if len(tbl_bytes) > 0:
            try:
                tbl_files[tbl_name] = tbl_bytes.decode("utf-8")
            except UnicodeDecodeError:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail=f"Restraints file '{tbl_name}' must be valid UTF-8 text.",
                )

    job_id = str(uuid.uuid4())

    record = await modal_manager.submit_job(
        job_id=job_id,
        pdb_files=pdb_files,
        tbl_files=tbl_files,
        sampling=sampling,
        refinement=refinement,
        gpu_type=gpu_type.value,
        gpu_device=gpu_device,
        gpu_platform=gpu_platform,
        job_name=job_name,
    )

    return JobSubmissionResponse(
        job_id=job_id,
        status=record.status,
        gpu_type=gpu_type.value,
        modal_call_id=record.modal_call_id,
        status_url=f"/api/v1/docking/{job_id}/status",
        results_url=f"/api/v1/results/{job_id}",
        created_at=record.created_at,
    )


@router.get(
    "/{job_id}/status",
    response_model=JobStatusResponse,
    summary="Get execution status and stage progress of a docking job",
)
async def get_docking_job_status(job_id: str) -> JobStatusResponse:
    """Retrieve the real-time execution stage and lifecycle status of a submitted job."""
    status_response = await modal_manager.get_status(job_id)
    if status_response is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job '{job_id}' not found in registry.",
        )
    return status_response


@router.post(
    "/{job_id}/cancel",
    summary="Cancel a running or queued docking job",
)
async def cancel_docking_job(job_id: str) -> dict[str, str]:
    """Terminate the remote Modal cloud GPU execution for the specified job."""
    success = await modal_manager.cancel_job(job_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job '{job_id}' not found in registry.",
        )
    return {
        "job_id": job_id,
        "status": JobStatus.CANCELLED.value,
        "message": "Docking job successfully cancelled on Modal compute plane.",
    }
