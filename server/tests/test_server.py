"""Integration tests for HADDOCK3 FastAPI Cloud Server."""

import io
from unittest.mock import MagicMock, patch

import pytest
from httpx import ASGITransport, AsyncClient

from server.client import modal_manager
from server.main import app

DUMMY_PDB_1 = (
    "ATOM      1  N   MET A   1      27.240  24.414  25.914  1.00 11.22           N\n"
    "ATOM      2  CA  MET A   1      26.505  25.267  26.857  1.00 12.01           C\n"
    "ATOM      3  N   GLY B   1      20.000  20.000  20.000  1.00 10.00           N\n"
    "END\n"
)

DUMMY_PDB_2 = (
    "ATOM      1  N   ALA B   1      15.120  18.234  12.441  1.00 10.50           N\n"
    "ATOM      2  CA  ALA B   1      14.300  19.100  13.200  1.00 10.90           C\n"
    "END\n"
)


@pytest.fixture(autouse=True)
def clean_registry():
    """Clear in-memory jobs before each test."""
    modal_manager.jobs.clear()


@pytest.mark.asyncio
async def test_root_endpoint():
    """Verify root index returns valid service metadata and documentation links."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        response = await client.get("/")
        assert response.status_code == 200
        data = response.json()
        assert "HADDOCK3 GPU" in data["title"]
        assert data["docs_url"] == "/docs"
        assert data["version"] == "3.0.0-gpu"


@pytest.mark.asyncio
async def test_system_health():
    """Verify Cloud Run health probe returns healthy status."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        response = await client.get("/api/v1/system/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert "timestamp" in data


@pytest.mark.asyncio
async def test_system_gpu():
    """Verify GPU metadata endpoint lists Modal accelerators."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        response = await client.get("/api/v1/system/gpu")
        assert response.status_code == 200
        data = response.json()
        assert "A100" in data["supported_gpus"]
        assert data["default_gpu"] == "A100"
        assert data["cloud_run_environment"] is True


@pytest.mark.asyncio
async def test_submit_validation_missing_molecule():
    """Submit must fail with 422 if either required molecule is missing."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        files = {
            "mol1_file": ("mol1.pdb", io.BytesIO(DUMMY_PDB_1.encode("utf-8")), "text/plain"),
        }
        response = await client.post("/api/v1/docking/submit", files=files)
        assert response.status_code == 422  # Missing mol2_file


@pytest.mark.asyncio
async def test_submit_validation_invalid_file_extension():
    """Submit must reject unsupported file extensions."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        files = {
            "mol1_file": ("mol1.pdb", io.BytesIO(DUMMY_PDB_1.encode("utf-8")), "text/plain"),
            "mol2_file": ("mol2.txt", io.BytesIO(DUMMY_PDB_2.encode("utf-8")), "text/plain"),
        }
        response = await client.post("/api/v1/docking/submit", files=files)
        assert response.status_code == 400
        assert "Unsupported file format" in response.json()["detail"]


@pytest.mark.asyncio
async def test_submit_validation_chain_not_found():
    """Submit must fail with 400 if user selects a chain that does not exist in PDB."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        files = {
            "mol1_file": ("mol1.pdb", io.BytesIO(DUMMY_PDB_1.encode("utf-8")), "text/plain"),
            "mol2_file": ("mol2.pdb", io.BytesIO(DUMMY_PDB_2.encode("utf-8")), "text/plain"),
        }
        data = {
            "mol1_chain": "Z",  # Chain Z doesn't exist in DUMMY_PDB_1
        }
        response = await client.post("/api/v1/docking/submit", files=files, data=data)
        assert response.status_code == 400
        assert "Specified chain 'Z' not found" in response.json()["detail"]


@pytest.mark.asyncio
async def test_submit_and_poll_workflow_with_params():
    """Test full submit -> status -> cancel lifecycle with separate molecule fields & params."""
    mock_call = MagicMock()
    mock_call.object_id = "mock-modal-call-777"
    mock_call.get.side_effect = TimeoutError()
    mock_call.get.aio = MagicMock(side_effect=TimeoutError())

    mock_fn = MagicMock()
    mock_fn.spawn.return_value = mock_call

    async def _async_spawn(**kwargs):
        _async_spawn.called_kwargs = kwargs
        return mock_call

    mock_fn.spawn.aio = _async_spawn

    with patch("server.client.modal.Function.from_name", return_value=mock_fn):
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as client:
            files = {
                "mol1_file": ("mol1.pdb", io.BytesIO(DUMMY_PDB_1.encode("utf-8")), "text/plain"),
                "mol2_file": ("mol2.pdb", io.BytesIO(DUMMY_PDB_2.encode("utf-8")), "text/plain"),
            }
            data = {
                "mol1_chain": "A",
                "mol1_kind": "Protein or Protein-Ligand",
                "mol1_coarse_grain": "false",
                "mol1_cyclic_peptide": "true",
                "mol2_chain": "B",
                "mol2_kind": "Protein or Protein-Ligand",
                "mol2_coarse_grain": "false",
                "mol2_cyclic_peptide": "false",
                "sampling": 10,
                "refinement": 5,
                "gpu_type": "A100",
                "job_name": "Portal Style Test",
            }
            # 1. Submit
            sub_resp = await client.post("/api/v1/docking/submit", files=files, data=data)
            assert sub_resp.status_code == 202
            sub_data = sub_resp.json()
            job_id = sub_data["job_id"]
            assert sub_data["status"] == "running"
            assert sub_data["gpu_type"] == "A100"
            assert sub_data["modal_call_id"] == "mock-modal-call-777"

            # Check that spawn received filtered mol1 PDB (only Chain A, not Chain B)
            called_kwargs = getattr(_async_spawn, "called_kwargs", {})
            saved_mol1_text = called_kwargs["pdb_files"]["mol1.pdb"]
            assert "MET A   1" in saved_mol1_text
            assert "GLY B   1" not in saved_mol1_text  # Chain B stripped out!
            # Check cyclic peptide param passed
            assert called_kwargs["mol_params"]["mol1"]["cyclic_peptide"] is True

            # 2. Poll Status (Running)
            stat_resp = await client.get(f"/api/v1/docking/{job_id}/status")
            assert stat_resp.status_code == 200
            stat_data = stat_resp.json()
            assert stat_data["job_id"] == job_id
            assert stat_data["status"] == "running"

            # 3. Check Results early (should return 409 Conflict)
            res_early = await client.get(f"/api/v1/results/{job_id}")
            assert res_early.status_code == 409
            assert "is currently running" in res_early.json()["detail"]

            # 4. Cancel
            cancel_resp = await client.post(f"/api/v1/docking/{job_id}/cancel")
            assert cancel_resp.status_code == 200
            assert cancel_resp.json()["status"] == "cancelled"
            mock_call.cancel.assert_called_once()


@pytest.mark.asyncio
async def test_completed_results_and_download():
    """Verify results parsing and PDB artifact file download for a completed job."""
    mock_results = {
        "status": "completed",
        "clusters": [
            {
                "cluster_rank": 1,
                "cluster_id": 1,
                "model_count": 4,
                "haddock_score": -89.42,
                "fnat": 0.82,
                "irmsd": 1.25,
                "lrmsd": 2.10,
                "dockq": 0.81,
            }
        ],
        "best_models": ["cluster1_1.pdb"],
        "pdb_artifacts": {
            "cluster1_1.pdb": DUMMY_PDB_1,
        },
    }

    mock_call = MagicMock()
    mock_call.object_id = "mock-modal-call-888"
    mock_call.get.return_value = mock_results

    async def _async_get(timeout=0):
        return mock_results

    mock_call.get.aio = _async_get

    mock_fn = MagicMock()
    mock_fn.spawn.return_value = mock_call

    async def _async_spawn(**kwargs):
        return mock_call

    mock_fn.spawn.aio = _async_spawn

    with patch("server.client.modal.Function.from_name", return_value=mock_fn):
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as client:
            files = {
                "mol1_file": ("mol1.pdb", io.BytesIO(DUMMY_PDB_1.encode("utf-8")), "text/plain"),
                "mol2_file": ("mol2.pdb", io.BytesIO(DUMMY_PDB_2.encode("utf-8")), "text/plain"),
            }
            # Submit
            sub_resp = await client.post("/api/v1/docking/submit", files=files)
            job_id = sub_resp.json()["job_id"]

            # Poll status (triggers completion and artifact saving)
            stat_resp = await client.get(f"/api/v1/docking/{job_id}/status")
            assert stat_resp.status_code == 200
            assert stat_resp.json()["status"] == "completed"

            # Fetch results
            res_resp = await client.get(f"/api/v1/results/{job_id}")
            assert res_resp.status_code == 200
            res_data = res_resp.json()
            assert len(res_data["clusters"]) == 1
            assert res_data["clusters"][0]["haddock_score"] == -89.42
            assert res_data["clusters"][0]["dockq"] == 0.81
            assert "cluster1_1.pdb" in res_data["best_models"]
            assert "cluster1_1.pdb" in res_data["download_urls"]

            # Download structure artifact
            dl_resp = await client.get(f"/api/v1/results/{job_id}/download/cluster1_1.pdb")
            assert dl_resp.status_code == 200
            assert dl_resp.text == DUMMY_PDB_1


@pytest.mark.asyncio
async def test_download_security_traversal_prevention():
    """Verify directory traversal is strictly blocked on artifact download."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        dl_resp = await client.get("/api/v1/results/some-job/download/..%2F..%2Fetc%2Fpasswd")
        assert dl_resp.status_code in (400, 404)
