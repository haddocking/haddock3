"""GPU acceleration utilities and device management for HADDOCK3."""

import os
import shutil
import subprocess
from pathlib import Path

from haddock import log


def is_torch_cuda_available() -> bool:
    """Check if PyTorch with CUDA support is available."""
    try:
        import torch

        return bool(torch.cuda.is_available())
    except (ImportError, RuntimeError, AttributeError):
        return False


def is_torch_mps_available() -> bool:
    """Check if PyTorch with Apple Silicon MPS support is available."""
    try:
        import torch

        return bool(
            hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
        )
    except (ImportError, RuntimeError, AttributeError):
        return False


def is_openmm_gpu_available() -> bool:
    """Check if OpenMM has CUDA or OpenCL platforms available."""
    try:
        import openmm

        platforms = [
            openmm.Platform.getPlatform(i).getName()
            for i in range(openmm.Platform.getNumPlatforms())
        ]
        return any(p in ("CUDA", "OpenCL") for p in platforms)
    except (ImportError, RuntimeError, AttributeError):
        return False


def is_gpu_available() -> bool:
    """Check whether any supported GPU acceleration platform is available.

    Returns:
        bool: True if CUDA, OpenCL, or MPS is detected and usable.
    """
    return (
        is_torch_cuda_available()
        or is_torch_mps_available()
        or is_openmm_gpu_available()
    )


def get_available_gpus() -> list[int]:
    """Retrieve list of available CUDA GPU device indices.

    Returns:
        list[int]: List of GPU indices, e.g. [0, 1]. Returns [0] if MPS is
            available on macOS, or an empty list if no GPU is found.
    """
    if is_torch_cuda_available():
        import torch

        return list(range(torch.cuda.device_count()))

    # Fallback to nvidia-smi if torch is not installed
    nvidia_smi = shutil.which("nvidia-smi")
    if nvidia_smi:
        try:
            out = subprocess.check_output(
                [nvidia_smi, "--query-gpu=index", "--format=csv,noheader"],
                text=True,
            )
            return [
                int(line.strip())
                for line in out.strip().splitlines()
                if line.strip().isdigit()
            ]
        except (subprocess.SubprocessError, OSError) as err:
            log.debug(f"nvidia-smi query failed: {err}")

    if is_torch_mps_available():
        return [0]

    return []


def get_gpu_device_name(device_id: int = 0) -> str:
    """Retrieve the model name for a specific GPU device.

    Args:
        device_id: GPU device index (default 0).

    Returns:
        str: Human-readable device name or 'CPU / Unknown'.
    """
    if is_torch_cuda_available():
        import torch

        try:
            return str(torch.cuda.get_device_name(device_id))
        except (RuntimeError, ValueError) as err:
            log.debug(f"Could not retrieve device name for GPU {device_id}: {err}")

    if is_torch_mps_available():
        return "Apple Silicon (MPS)"

    return "CPU / Unknown"


def resolve_gpu_platform(requested_platform: str = "auto") -> str:
    """Resolve acceleration backend platform based on availability and preference.

    Args:
        requested_platform: Requested platform choice: 'auto', 'cuda', 'opencl',
            'mps', or 'cpu'.

    Returns:
        str: Selected platform string ('cuda', 'opencl', 'mps', or 'cpu').
    """
    platform_req = requested_platform.lower()

    if platform_req == "cpu":
        return "cpu"

    if platform_req == "cuda":
        if is_torch_cuda_available() or is_openmm_gpu_available():
            return "cuda"
        log.warning("CUDA platform requested but unavailable. Falling back to CPU.")
        return "cpu"

    if platform_req == "mps":
        if is_torch_mps_available():
            return "mps"
        log.warning("MPS platform requested but unavailable. Falling back to CPU.")
        return "cpu"

    if platform_req == "opencl":
        if is_openmm_gpu_available():
            return "opencl"
        log.warning("OpenCL platform requested but unavailable. Falling back to CPU.")
        return "cpu"

    # 'auto' mode: Prioritize CUDA > MPS > OpenCL > CPU
    if is_torch_cuda_available():
        return "cuda"
    if is_openmm_gpu_available():
        return "cuda"
    if is_torch_mps_available():
        return "mps"

    return "cpu"


def get_best_available_device(
    preferred_platform: str = "auto",
    gpu_devices: list[int] | None = None,
) -> str:
    """Return the torch/OpenMM compatible device string based on available hardware.

    Args:
        preferred_platform: Preferred platform ('auto', 'cuda', 'mps', 'cpu').
        gpu_devices: Optional list of specific GPU device IDs to select from.

    Returns:
        str: Device identifier string ('cuda', 'cuda:0', 'mps', or 'cpu').
    """
    platform = resolve_gpu_platform(preferred_platform)
    if platform == "cuda":
        if gpu_devices:
            return f"cuda:{gpu_devices[0]}"
        return "cuda"
    if platform == "mps":
        return "mps"
    return "cpu"



class GPUDevicePool:
    """Manages assignment of GPU devices to parallel workers."""

    def __init__(self, devices: list[int] | None = None) -> None:
        """Initialize the GPU device pool.

        Args:
            devices: List of GPU device indices to utilize. If None or empty,
                auto-detects available devices or defaults to [0].
        """
        if devices:
            self.devices: list[int] = list(devices)
        else:
            avail = get_available_gpus()
            self.devices = avail if avail else [0]

    def get_device(self, worker_index: int) -> int:
        """Assign a GPU device index for a given worker via round-robin.

        Args:
            worker_index: Integer index of the worker process.

        Returns:
            int: Assigned GPU device index.
        """
        if not self.devices:
            return 0
        return self.devices[worker_index % len(self.devices)]

    @staticmethod
    def set_cuda_visible_device(device_id: int | None) -> None:
        """Set CUDA_VISIBLE_DEVICES for the current process.

        Args:
            device_id: GPU device index to expose, or None to leave unset.
        """
        if device_id is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = str(device_id)
            log.debug(
                f"Assigned CUDA_VISIBLE_DEVICES={device_id} to process {os.getpid()}"
            )


def is_nvidia_mps_running() -> bool:
    """Check if the NVIDIA Multi-Process Service (MPS) control daemon is active.

    Returns:
        bool: True if nvidia-cuda-mps-control is running.
    """
    mps_pipe_dir = os.environ.get("CUDA_MPS_PIPE_DIRECTORY", "/tmp/nvidia-mps")
    control_socket = Path(mps_pipe_dir, "control")
    return control_socket.exists()


def start_nvidia_mps_daemon() -> bool:
    """Start NVIDIA CUDA MPS control daemon for multiplexing workers onto a GPU.

    Returns:
        bool: True if daemon started or was already active, False on failure.
    """
    if is_nvidia_mps_running():
        return True

    mps_exec = shutil.which("nvidia-cuda-mps-control")
    if not mps_exec:
        log.debug("nvidia-cuda-mps-control not found on system PATH.")
        return False

    try:
        subprocess.run(
            [mps_exec, "-d"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        log.info("Started NVIDIA CUDA Multi-Process Service (MPS) daemon.")
        return True
    except (subprocess.SubprocessError, OSError) as err:
        log.warning(f"Could not start NVIDIA CUDA MPS daemon: {err}")
        return False


def stop_nvidia_mps_daemon() -> bool:
    """Shut down NVIDIA CUDA MPS control daemon.

    Returns:
        bool: True if stopped or not running, False on failure.
    """
    if not is_nvidia_mps_running():
        return True

    mps_exec = shutil.which("nvidia-cuda-mps-control")
    if not mps_exec:
        return True

    try:
        subprocess.run(
            [mps_exec],
            input=b"quit\n",
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        log.info("Stopped NVIDIA CUDA Multi-Process Service (MPS) daemon.")
        return True
    except (subprocess.SubprocessError, OSError) as err:
        log.warning(f"Could not stop NVIDIA CUDA MPS daemon: {err}")
        return False

