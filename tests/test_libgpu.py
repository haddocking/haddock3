"""Unit tests for haddock.libs.libgpu."""

import os
from unittest.mock import MagicMock, patch

from haddock.libs.libgpu import (
    GPUDevicePool,
    get_available_gpus,
    get_gpu_device_name,
    is_gpu_available,
    is_openmm_gpu_available,
    is_torch_cuda_available,
    is_torch_mps_available,
    resolve_gpu_platform,
)


def test_is_torch_cuda_available_true():
    """Test torch CUDA detection when available."""
    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = True
    with patch.dict("sys.modules", {"torch": mock_torch}):
        assert is_torch_cuda_available() is True


def test_is_torch_cuda_available_false():
    """Test torch CUDA detection when unavailable."""
    with patch.dict("sys.modules", {"torch": None}):
        assert is_torch_cuda_available() is False


def test_is_torch_mps_available():
    """Test torch MPS detection."""
    mock_torch = MagicMock()
    mock_torch.backends.mps.is_available.return_value = True
    with patch.dict("sys.modules", {"torch": mock_torch}):
        assert is_torch_mps_available() is True


def test_is_openmm_gpu_available():
    """Test OpenMM platform detection for CUDA and OpenCL."""
    mock_openmm = MagicMock()
    mock_platform_cuda = MagicMock()
    mock_platform_cuda.getName.return_value = "CUDA"
    mock_openmm.Platform.getNumPlatforms.return_value = 1
    mock_openmm.Platform.getPlatform.return_value = mock_platform_cuda
    with patch.dict("sys.modules", {"openmm": mock_openmm}):
        assert is_openmm_gpu_available() is True


def test_is_gpu_available_fallback():
    """Test is_gpu_available returns boolean without error even if no packages installed."""
    result = is_gpu_available()
    assert isinstance(result, bool)


def test_get_available_gpus():
    """Test device enumeration."""
    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = True
    mock_torch.cuda.device_count.return_value = 2
    with patch.dict("sys.modules", {"torch": mock_torch}):
        devices = get_available_gpus()
        assert devices == [0, 1]


def test_get_gpu_device_name():
    """Test retrieving GPU device name."""
    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = True
    mock_torch.cuda.get_device_name.return_value = "NVIDIA A100-SXM4-80GB"
    with patch.dict("sys.modules", {"torch": mock_torch}):
        name = get_gpu_device_name(0)
        assert "A100" in name


def test_resolve_gpu_platform():
    """Test resolving platform strings."""
    assert resolve_gpu_platform("cpu") == "cpu"

    # Test auto fallback to cpu when no gpu
    with (
        patch("haddock.libs.libgpu.is_torch_cuda_available", return_value=False),
        patch("haddock.libs.libgpu.is_openmm_gpu_available", return_value=False),
        patch("haddock.libs.libgpu.is_torch_mps_available", return_value=False),
    ):
        assert resolve_gpu_platform("auto") == "cpu"
        assert resolve_gpu_platform("cuda") == "cpu"


def test_gpu_device_pool():
    """Test GPUDevicePool round-robin device distribution."""
    pool = GPUDevicePool([0, 1, 2])
    assert pool.get_device(0) == 0
    assert pool.get_device(1) == 1
    assert pool.get_device(2) == 2
    assert pool.get_device(3) == 0
    assert pool.get_device(4) == 1

    orig_cuda = os.environ.get("CUDA_VISIBLE_DEVICES")
    try:
        GPUDevicePool.set_cuda_visible_device(1)
        assert os.environ.get("CUDA_VISIBLE_DEVICES") == "1"
    finally:
        if orig_cuda is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = orig_cuda
        else:
            os.environ.pop("CUDA_VISIBLE_DEVICES", None)


def test_get_best_available_device():
    """Test get_best_available_device under different hardware availability states."""
    from haddock.libs.libgpu import get_best_available_device

    assert get_best_available_device("cpu") == "cpu"
    with patch("haddock.libs.libgpu.is_torch_cuda_available", return_value=True):
        assert get_best_available_device("auto") == "cuda"
        assert get_best_available_device("cuda", gpu_devices=[2]) == "cuda:2"


def test_nvidia_mps_daemon_handling():
    """Test NVIDIA MPS daemon start and stop functions."""
    from haddock.libs.libgpu import (
        is_nvidia_mps_running,
        start_nvidia_mps_daemon,
        stop_nvidia_mps_daemon,
    )

    assert isinstance(is_nvidia_mps_running(), bool)

    # When executable not present
    with patch("shutil.which", return_value=None):
        assert start_nvidia_mps_daemon() is False
        assert stop_nvidia_mps_daemon() is True

    # When executable present and successful
    with (
        patch("shutil.which", return_value="/usr/bin/nvidia-cuda-mps-control"),
        patch("subprocess.run") as mock_run,
    ):
        mock_run.return_value = MagicMock(returncode=0)
        with patch("haddock.libs.libgpu.is_nvidia_mps_running", return_value=False):
            assert start_nvidia_mps_daemon() is True
            assert mock_run.called

