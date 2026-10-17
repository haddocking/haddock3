"""Tests for accelerated FCC computation in libfcc_gpu."""

import math
import tempfile
from pathlib import Path

import numpy as np

from haddock.libs.libfcc import calculate_pairwise_matrix
from haddock.libs.libfcc_gpu import (
    calculate_pairwise_matrix_gpu,
    write_fcc_matrix_file,
)


def test_empty_and_single_model():
    """Test FCC calculation on edge cases (empty and single model)."""
    assert calculate_pairwise_matrix_gpu([]) == []
    assert calculate_pairwise_matrix_gpu([set()]) == []
    assert calculate_pairwise_matrix_gpu([{1, 2, 3}]) == []


def test_zero_contacts():
    """Test FCC calculation when models have no contacts."""
    contacts = [set(), set(), set()]
    res = calculate_pairwise_matrix_gpu(contacts)
    assert len(res) == 3
    assert res == [(1, 2, 0.0, 0.0), (1, 3, 0.0, 0.0), (2, 3, 0.0, 0.0)]


def test_identical_contacts():
    """Test FCC calculation when models have identical contacts."""
    contacts = [{10, 20, 30}, {10, 20, 30}]
    res = calculate_pairwise_matrix_gpu(contacts)
    assert len(res) == 1
    assert res[0] == (1, 2, 1.0, 1.0)


def test_disjoint_contacts():
    """Test FCC calculation when models share no contacts."""
    contacts = [{1, 2}, {3, 4}, {5, 6}]
    res = calculate_pairwise_matrix_gpu(contacts)
    assert len(res) == 3
    for entry in res:
        assert entry[2] == 0.0
        assert entry[3] == 0.0


def test_equivalence_with_cpu_loop():
    """Verify numerical equivalence between CPU loop and accelerated matrix."""
    rng = np.random.default_rng(12345)
    pool = list(range(500))
    contacts = [
        set(rng.choice(pool, size=rng.integers(20, 80), replace=False))
        for _ in range(50)
    ]

    cpu_ref = list(calculate_pairwise_matrix(contacts, ignore_chain=False))
    fast_res = calculate_pairwise_matrix_gpu(contacts, chunk_size=16)

    assert len(cpu_ref) == len(fast_res)
    for (i1, k1, f1, f2), (i2, k2, g1, g2) in zip(cpu_ref, fast_res):
        assert i1 == i2
        assert k1 == k2
        assert math.isclose(f1, g1, abs_tol=1e-5)
        assert math.isclose(f2, g2, abs_tol=1e-5)


def test_write_fcc_matrix_file():
    """Verify write_fcc_matrix_file generates exact expected file format."""
    matrix = [
        (1, 2, 0.5, 0.666666),
        (1, 3, 0.25, 0.25),
        (2, 3, 0.333333, 0.25),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        out_f = Path(tmpdir, "test.matrix")
        write_fcc_matrix_file(matrix, out_f, buffer_lines=2)
        assert out_f.exists()

        content = out_f.read_text().strip().splitlines()
        assert len(content) == 3
        assert content[0] == "1 2 0.50 0.667"
        assert content[1] == "1 3 0.25 0.250"
        assert content[2] == "2 3 0.33 0.250"
