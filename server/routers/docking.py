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


def filter_pdb_by_chain(pdb_text: str, chain: str) -> str:
    """Filter PDB file contents to keep only atoms belonging to the specified chain.

    If chain is 'All', '*' or empty, the content is returned unmodified.
    """
    clean_chain = chain.strip()
    if clean_chain.lower() in ("all", "*", ""):
        return pdb_text

    target_chain = clean_chain.upper()
    filtered_lines = []
    found_target_atoms = False

    for line in pdb_text.splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            chain_id = line[21].strip() if len(line) > 21 else ""
            if chain_id.upper() == target_chain:
                filtered_lines.append(line)
                found_target_atoms = True
        elif line.startswith(("TER", "ANISOU")):
            chain_id = line[21].strip() if len(line) > 21 else ""
            if not chain_id or chain_id.upper() == target_chain:
                filtered_lines.append(line)
        else:
            # Preserve headers, SEQRES, CRYST1, CONECT, END
            filtered_lines.append(line)

    if not found_target_atoms:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Specified chain '{chain}' not found in molecular structure.",
        )

    return "".join(filtered_lines)


async def _read_and_validate_structure(
    file: UploadFile,
    chain: str,
    default_name: str,
) -> tuple[str, str]:
    """Validate structure file format, read bytes, and apply chain filter."""
    fname = file.filename or default_name
    valid_exts = (".pdb", ".ent", ".cif")
    if not any(fname.lower().endswith(ext) for ext in valid_exts):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Unsupported file format for '{fname}'. Must be .pdb, .ent, or .cif.",
        )

    max_bytes = settings.MAX_UPLOAD_SIZE_MB * 1024 * 1024
    content_bytes = await file.read()
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
        raw_text = content_bytes.decode("utf-8")
    except UnicodeDecodeError:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"File '{fname}' must be valid UTF-8 text.",
        )

    filtered_text = filter_pdb_by_chain(raw_text, chain)
    return fname, filtered_text


@router.post(
    "/submit",
    response_model=JobSubmissionResponse,
    status_code=status.HTTP_202_ACCEPTED,
    summary="Submit a molecular docking calculation to Modal GPU",
    description=(
        "Submits 2 molecular structures with dedicated per-molecule configuration "
        "(chain selection, molecule kind, coarse-graining, cyclic peptide detection) "
        "and optional restraint files for accelerated docking on NVIDIA A100 GPUs via Modal."
    ),
)
async def submit_docking_job(
    # Molecule 1 - input
    mol1_file: UploadFile = File(
        ...,
        description="PDB or mmCIF structure to submit for Molecule 1.",
    ),
    mol1_chain: str = Form(
        "All",
        description="Which chain of Molecule 1 must be used? ('All' or specific chain ID e.g. 'A').",
    ),
    mol1_kind: str = Form(
        "Protein or Protein-Ligand",
        description="What kind of molecule are you docking? (e.g. 'Protein or Protein-Ligand', 'DNA', 'RNA', 'Small Molecule').",
    ),
    mol1_coarse_grain: bool = Form(
        False,
        description="Do you want to coarse-grain Molecule 1? Convert all-atom structure into Martini coarse-grained.",
    ),
    mol1_cyclic_peptide: bool = Form(
        False,
        description="Is Molecule 1 a cyclic peptide? HADDOCK will generate a peptide bond between N- and C-termini.",
    ),
    # Molecule 2 - input
    mol2_file: UploadFile = File(
        ...,
        description="PDB or mmCIF structure to submit for Molecule 2.",
    ),
    mol2_chain: str = Form(
        "All",
        description="Which chain of Molecule 2 must be used? ('All' or specific chain ID e.g. 'B').",
    ),
    mol2_kind: str = Form(
        "Protein or Protein-Ligand",
        description="What kind of molecule are you docking? (e.g. 'Protein or Protein-Ligand', 'DNA', 'RNA', 'Small Molecule').",
    ),
    mol2_coarse_grain: bool = Form(
        False,
        description="Do you want to coarse-grain Molecule 2? Convert all-atom structure into Martini coarse-grained.",
    ),
    mol2_cyclic_peptide: bool = Form(
        False,
        description="Is Molecule 2 a cyclic peptide? HADDOCK will generate a peptide bond between N- and C-termini.",
    ),
    # Optional Molecule 3 (for multi-body complexes)
    mol3_file: Optional[UploadFile] = File(
        None,
        description="Optional PDB or mmCIF structure for Molecule 3 (multi-body docking).",
    ),
    mol3_chain: Optional[str] = Form(
        "All",
        description="Which chain of Molecule 3 must be used?",
    ),
    # Restraints and workflow configuration
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
    pdb_files: dict[str, str] = {}

    # Process Molecule 1
    m1_name, m1_text = await _read_and_validate_structure(
        mol1_file, chain=mol1_chain, default_name="mol1.pdb"
    )
    pdb_files[m1_name] = m1_text

    # Process Molecule 2
    m2_name, m2_text = await _read_and_validate_structure(
        mol2_file, chain=mol2_chain, default_name="mol2.pdb"
    )
    # Ensure distinct filename in case both uploaded as 'model.pdb'
    if m2_name == m1_name:
        m2_name = f"mol2_{m2_name}"
    pdb_files[m2_name] = m2_text

    # Process Molecule 3 if provided
    if mol3_file is not None and bool(mol3_file.filename):
        m3_name, m3_text = await _read_and_validate_structure(
            mol3_file, chain=mol3_chain or "All", default_name="mol3.pdb"
        )
        if m3_name in pdb_files:
            m3_name = f"mol3_{m3_name}"
        pdb_files[m3_name] = m3_text

    # Process Restraints
    tbl_files: dict[str, str] = {}
    if restraints is not None and bool(restraints.filename):
        max_bytes = settings.MAX_UPLOAD_SIZE_MB * 1024 * 1024
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

    mol_params = {
        "mol1": {
            "chain": mol1_chain,
            "kind": mol1_kind,
            "coarse_grain": mol1_coarse_grain,
            "cyclic_peptide": mol1_cyclic_peptide,
        },
        "mol2": {
            "chain": mol2_chain,
            "kind": mol2_kind,
            "coarse_grain": mol2_coarse_grain,
            "cyclic_peptide": mol2_cyclic_peptide,
        },
    }

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
        mol_params=mol_params,
    )

    return JobSubmissionResponse(
        job_id=job_id,
        status=record.status,
        gpu_type=gpu_type.value,
        modal_call_id=record.modal_call_id,
        status_url=f"/api/v1/docking/{job_id}/status",
        results_url=f"/api/v1/results/{job_id}",
        created_at=record.created_at,
        error_message=record.error_message,
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
