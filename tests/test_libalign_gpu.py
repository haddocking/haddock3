"""Unit tests for GPU-accelerated and tensorized RMSD algorithms."""

import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from haddock.libs.libalign_gpu import (
    compute_rmsd_matrix,
    compute_rmsd_matrix_numpy,
    write_rmsd_matrix_file,
)


@pytest.fixture
def sample_structures():
    """Create sample coordinates for 3 models with 10 atoms."""
    p = np.array(
        [
            [-3.811, -0.120, 7.200],
            [-5.913, 8.536, 6.286],
            [12.839, 8.507, 4.587],
            [4.856, -3.766, -10.577],
            [-1.922, -4.171, -6.596],
            [-12.186, -4.807, -10.777],
            [0.308, -7.500, 8.632],
            [-0.080, -7.218, 2.967],
            [9.232, 18.419, 8.949],
            [-1.748, -8.395, 7.604],
        ],
        dtype=np.float32,
    )

    q = np.array(
        [
            [6.095, -5.630, 8.033],
            [0.508, 1.733, 10.833],
            [7.427, 7.670, -5.757],
            [-12.593, -6.526, 5.217],
            [-6.309, -4.756, 2.355],
            [-3.664, -8.990, -8.014],
            [11.789, -10.623, 3.966],
            [6.895, -12.338, 1.268],
            [5.764, 16.110, 1.344],
            [10.306, -12.226, 5.028],
        ],
        dtype=np.float32,
    )

    # Third structure is shifted version of P
    r = p + np.array([2.0, -1.0, 3.0], dtype=np.float32)

    coords = np.stack([p, q, r])
    return coords


def test_compute_rmsd_matrix_numpy(sample_structures):
    """Test NumPy batched RMSD calculation."""
    i_idx, j_idx, rmsds = compute_rmsd_matrix_numpy(sample_structures)
    assert len(i_idx) == 3
    assert len(j_idx) == 3
    assert len(rmsds) == 3

    # Pairs: (0, 1), (0, 2), (1, 2)
    assert (i_idx[0], j_idx[0]) == (0, 1)
    assert (i_idx[1], j_idx[1]) == (0, 2)
    assert (i_idx[2], j_idx[2]) == (1, 2)

    # Model 0 and Model 2 are identical structures shifted in space -> RMSD must be 0
    assert rmsds[1] == pytest.approx(0.0, abs=1e-5)
    # Model 0 and Model 1 pairwise Kabsch RMSD
    assert rmsds[0] == pytest.approx(7.591, abs=1e-2)


def test_write_rmsd_matrix_file(sample_structures):
    """Test writing RMSD matrix in HADDOCK3 format."""
    i_idx, j_idx, rmsds = compute_rmsd_matrix_numpy(sample_structures)
    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = Path(tmpdir, "rmsd.matrix")
        write_rmsd_matrix_file(out_path, i_idx, j_idx, rmsds)
        assert out_path.exists()

        lines = out_path.read_text().splitlines()
        assert len(lines) == 3
        # Indices in output file must be 1-based
        assert lines[0].startswith("1 2 ")
        assert lines[1].startswith("1 3 0.000")
        assert lines[2].startswith("2 3 ")


def test_compute_rmsd_matrix_fallback(sample_structures):
    """Test compute_rmsd_matrix falls back to numpy when torch is not available."""
    with patch(
        "haddock.libs.libalign_gpu.compute_rmsd_matrix_torch", side_effect=ImportError
    ):
        _, _, rmsds = compute_rmsd_matrix(sample_structures, use_gpu=True)
        assert len(rmsds) == 3
        assert rmsds[1] == pytest.approx(0.0, abs=1e-5)
