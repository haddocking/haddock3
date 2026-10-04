"""Accelerated Fraction of Common Contacts (FCC) matrix computation.

This module provides batched, tensorized, and sparse matrix implementations of the
pairwise Fraction of Common Contacts (FCC) algorithm used in clustering. It
supports GPU acceleration via PyTorch (CUDA, MPS) and CPU acceleration via SciPy
sparse matrix operations, eliminating the O(N^2) Python loop bottleneck.
"""

from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

from haddock import log

try:
    import torch

    TORCH_AVAILABLE = True
except ImportError:
    torch = None  # type: ignore[assignment]
    TORCH_AVAILABLE = False

try:
    from scipy.sparse import csr_matrix

    SCIPY_AVAILABLE = True
except ImportError:
    csr_matrix = None  # type: ignore[assignment]
    SCIPY_AVAILABLE = False


def calculate_pairwise_matrix_gpu(
    contacts: Sequence[Iterable[int]],
    device: str = "cpu",
    chunk_size: int = 2048,
) -> list[tuple[int, int, float, float]]:
    """Calculates pairwise fraction of common contacts using GPU or sparse tensors.

    Given N contact sets across complexes, computes the asymmetric pairwise FCC
    matrix:
        fcc(i, k) = |contacts[i] ∩ contacts[k]| / |contacts[i]|
        fcc(k, i) = |contacts[i] ∩ contacts[k]| / |contacts[k]|

    The computation represents contact sets as a binary matrix A of shape
    (N, total_unique_contacts) and evaluates intersection sizes via matrix
    multiplication M = A * A^T using GPU Tensor Cores or compiled SciPy sparse
    routines.

    Args:
        contacts: Sequence of length N containing collections of contact
            identifiers (integers representing residue-residue contacts).
        device: Computational device ('cpu', 'cuda', 'cuda:0', or 'mps').
        chunk_size: Number of rows to process per chunk to bound VRAM usage.

    Returns:
        A list of tuples in the format (i, k, fcc_ik, fcc_ki), where i and k
        are 1-indexed (1 <= i < k <= N) in upper-triangular order.
    """
    num_models = len(contacts)
    if num_models < 2:
        return []

    # Map unique contact IDs to contiguous column indices
    unique_contacts: dict[int, int] = {}
    total_unique = 0
    rows: list[int] = []
    cols: list[int] = []

    for i, con in enumerate(contacts):
        for c in con:
            idx = unique_contacts.get(c)
            if idx is None:
                idx = total_unique
                unique_contacts[c] = idx
                total_unique += 1
            rows.append(i)
            cols.append(idx)

    # If no contacts exist in any structure, all pairwise FCC values are 0.0
    if total_unique == 0:
        empty_res: list[tuple[int, int, float, float]] = []
        for i in range(num_models):
            for k in range(i + 1, num_models):
                empty_res.append((i + 1, k + 1, 0.0, 0.0))
        return empty_res

    # Pre-calculate inverse contact counts for normalization
    lens = np.array(
        [
            len(con) if hasattr(con, "__len__") else len(list(con))
            for con in contacts
        ],
        dtype=np.float64,
    )
    inv_lens = np.zeros_like(lens)
    nonzero = lens > 0
    inv_lens[nonzero] = 1.0 / lens[nonzero]

    # Attempt PyTorch execution if requested and available
    use_torch = TORCH_AVAILABLE and (
        device.startswith("cuda") or device == "mps" or not SCIPY_AVAILABLE
    )
    if use_torch and torch is not None:
        try:
            return _calculate_fcc_torch(
                num_models=num_models,
                total_unique=total_unique,
                rows=rows,
                cols=cols,
                inv_lens=inv_lens,
                device=device,
                chunk_size=chunk_size,
            )
        except (RuntimeError, ValueError, TypeError) as err:
            log.warning(
                f"PyTorch FCC calculation failed ({err}); falling back to CPU sparse."
            )

    # Fallback to SciPy sparse matrix on CPU
    if SCIPY_AVAILABLE and csr_matrix is not None:
        return _calculate_fcc_scipy(
            num_models=num_models,
            total_unique=total_unique,
            rows=rows,
            cols=cols,
            inv_lens=inv_lens,
            chunk_size=chunk_size,
        )

    # Pure Python fallback
    return _calculate_fcc_python_fallback(contacts, inv_lens)


def _calculate_fcc_torch(
    num_models: int,
    total_unique: int,
    rows: list[int],
    cols: list[int],
    inv_lens: np.ndarray,
    device: str,
    chunk_size: int,
) -> list[tuple[int, int, float, float]]:
    """Internal PyTorch tensor implementation for FCC calculation."""
    dev = torch.device(
        device if (torch.cuda.is_available() or device == "mps") else "cpu"
    )
    inv_lens_t = torch.from_numpy(inv_lens).to(device=dev, dtype=torch.float32)

    # Build binary matrix A of shape (N, total_unique)
    A = torch.zeros((num_models, total_unique), dtype=torch.float32, device=dev)
    row_t = torch.tensor(rows, dtype=torch.long, device=dev)
    col_t = torch.tensor(cols, dtype=torch.long, device=dev)
    A[row_t, col_t] = 1.0

    results: list[tuple[int, int, float, float]] = []

    # Process in row chunks to bound VRAM usage
    for i_start in range(0, num_models, chunk_size):
        i_end = min(i_start + chunk_size, num_models)
        A_chunk = A[i_start:i_end]
        # M_chunk has shape (chunk_size, N)
        M_chunk = torch.matmul(A_chunk, A.T)

        for r in range(i_end - i_start):
            i = i_start + r
            if i + 1 >= num_models:
                continue
            cols_idx = slice(i + 1, num_models)
            cc = M_chunk[r, cols_idx]
            fcc_ik = (cc * inv_lens_t[i]).cpu().tolist()
            fcc_ki = (cc * inv_lens_t[cols_idx]).cpu().tolist()

            k_indices = range(i + 2, num_models + 1)
            for k_val, f1, f2 in zip(k_indices, fcc_ik, fcc_ki):
                results.append((i + 1, k_val, float(f1), float(f2)))

    return results


def _calculate_fcc_scipy(
    num_models: int,
    total_unique: int,
    rows: list[int],
    cols: list[int],
    inv_lens: np.ndarray,
    chunk_size: int,
) -> list[tuple[int, int, float, float]]:
    """Internal SciPy sparse CSR implementation for FCC calculation."""
    data = np.ones(len(rows), dtype=np.float32)
    A = csr_matrix(
        (data, (rows, cols)), shape=(num_models, total_unique), dtype=np.float32
    )

    results: list[tuple[int, int, float, float]] = []
    for i_start in range(0, num_models, chunk_size):
        i_end = min(i_start + chunk_size, num_models)
        A_chunk = A[i_start:i_end]
        M_chunk = A_chunk.dot(A.T).toarray()

        for r in range(i_end - i_start):
            i = i_start + r
            if i + 1 >= num_models:
                continue
            cols_idx = slice(i + 1, num_models)
            cc = M_chunk[r, cols_idx]
            fcc_ik = cc * inv_lens[i]
            fcc_ki = cc * inv_lens[cols_idx]

            k_indices = range(i + 2, num_models + 1)
            for k_val, f1, f2 in zip(k_indices, fcc_ik, fcc_ki):
                results.append((i + 1, k_val, float(f1), float(f2)))

    return results


def _calculate_fcc_python_fallback(
    contacts: Sequence[Iterable[int]],
    inv_lens: np.ndarray,
) -> list[tuple[int, int, float, float]]:
    """Pure Python fallback for FCC calculation."""
    sets = [set(c) for c in contacts]
    results: list[tuple[int, int, float, float]] = []
    num_models = len(sets)

    for i in range(num_models):
        for k in range(i + 1, num_models):
            cc = len(sets[i].intersection(sets[k]))
            fcc = float(cc * inv_lens[i])
            fcc_v = float(cc * inv_lens[k])
            results.append((i + 1, k + 1, fcc, fcc_v))

    return results


def write_fcc_matrix_file(
    matrix: list[tuple[int, int, float, float]],
    output_path: Path | str,
    buffer_lines: int = 50000,
) -> None:
    """Writes an FCC pairwise matrix to disk with high-performance buffering.

    Args:
        matrix: List of tuples (i, k, fcc_ik, fcc_ki).
        output_path: Path to the target fcc.matrix file.
        buffer_lines: Number of text lines to accumulate before flushing.
    """
    buffer: list[str] = []
    with open(output_path, "w") as fh:
        for row in matrix:
            buffer.append(f"{row[0]} {row[1]} {row[2]:.2f} {row[3]:.3f}\n")
            if len(buffer) >= buffer_lines:
                fh.write("".join(buffer))
                buffer.clear()
        if buffer:
            fh.write("".join(buffer))
