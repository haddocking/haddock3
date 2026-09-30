"""Unit tests for GPU acceleration in OpenMM refinement module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.refinement.openmm.openmm import OPENMM


@pytest.fixture
def mock_openmm_job():
    """Create a mock OPENMM instance for unit testing."""
    params = {
        "use_gpu": True,
        "gpu_platform": "cuda",
        "constraints": "HBonds",
    }
    directory_dict = {
        "pdbfixer": "pdbfixer",
        "solvation_boxes": "solvation_boxes",
        "intermediates": "intermediates",
        "md_raw_output": "md_raw_output",
        "openmm_output": "openmm_output",
        "simulation_stats": "simulation_stats",
    }
    model = PDBFile(Path("dummy.pdb"))
    with patch("haddock.modules.refinement.openmm.openmm.HBonds", "HBonds"):
        job = OPENMM(
            identificator=1,
            model=model,
            path=Path("."),
            directory_dict=directory_dict,
            params=params,
        )
    return job


def test_openmm_create_simulation_cuda(mock_openmm_job):
    """Test OPENMM._create_simulation creates simulation on CUDA platform."""
    mock_platform = MagicMock()
    mock_simulation = MagicMock()
    mock_openmm_module = MagicMock()
    mock_openmm_module.Platform.getPlatformByName.return_value = mock_platform

    mock_openmm_job.gpu_device = 1
    mock_openmm_job.params["use_gpu"] = True
    mock_openmm_job.params["gpu_platform"] = "cuda"

    with (
        patch.dict("sys.modules", {"openmm": mock_openmm_module}),
        patch(
            "haddock.modules.refinement.openmm.openmm.Simulation",
            return_value=mock_simulation,
        ) as mock_sim_class,
        patch("haddock.libs.libgpu.is_torch_cuda_available", return_value=True),
    ):
        topology = MagicMock()
        system = MagicMock()
        integrator = MagicMock()

        sim = mock_openmm_job._create_simulation(topology, system, integrator)

        assert sim == mock_simulation
        mock_openmm_module.Platform.getPlatformByName.assert_called_with("CUDA")
        mock_sim_class.assert_called_once_with(
            topology,
            system,
            integrator,
            platform=mock_platform,
            platformProperties={"Precision": "mixed", "DeviceIndex": "1"},
        )


def test_openmm_create_simulation_fallback_on_error(mock_openmm_job):
    """Test OPENMM._create_simulation falls back gracefully if platform fails."""
    mock_openmm_module = MagicMock()
    mock_openmm_module.Platform.getPlatformByName.side_effect = RuntimeError(
        "CUDA driver error"
    )
    mock_simulation = MagicMock()

    mock_openmm_job.params["use_gpu"] = True
    mock_openmm_job.params["gpu_platform"] = "cuda"

    with (
        patch.dict("sys.modules", {"openmm": mock_openmm_module}),
        patch(
            "haddock.modules.refinement.openmm.openmm.Simulation",
            return_value=mock_simulation,
        ) as mock_sim_class,
        patch("haddock.libs.libgpu.is_torch_cuda_available", return_value=True),
    ):
        topology = MagicMock()
        system = MagicMock()
        integrator = MagicMock()

        sim = mock_openmm_job._create_simulation(topology, system, integrator)

        assert sim == mock_simulation
        mock_sim_class.assert_called_once_with(topology, system, integrator)


def test_openmm_create_simulation_cpu_default(mock_openmm_job):
    """Test OPENMM._create_simulation uses standard platform when use_gpu is False."""
    mock_openmm_job.params["use_gpu"] = False
    mock_simulation = MagicMock()

    with patch(
        "haddock.modules.refinement.openmm.openmm.Simulation",
        return_value=mock_simulation,
    ) as mock_sim_class:
        topology = MagicMock()
        system = MagicMock()
        integrator = MagicMock()

        sim = mock_openmm_job._create_simulation(topology, system, integrator)

        assert sim == mock_simulation
        mock_sim_class.assert_called_once_with(topology, system, integrator)
