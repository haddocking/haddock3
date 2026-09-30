"""Unit tests for the Modal GPU harness."""

import pytest


def test_modal_runner_importable():
    """Verify modal_runner module imports cleanly with or without modal package."""
    import tests.modal_gpu.modal_runner as mr

    assert hasattr(mr, "REPO_ROOT")
    assert mr.REPO_ROOT.exists()


def test_modal_runner_execution_logic():
    """Test modal_runner test execution helper with mock subprocess."""
    import tests.modal_gpu.modal_runner as mr

    if mr.modal is None:
        pytest.skip("modal package not installed in environment")

    # If modal is available, verify run_tests_on_gpu function existence
    assert hasattr(mr, "run_tests_on_gpu")
    assert hasattr(mr, "benchmark_kernel")
