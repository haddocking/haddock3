"""Tensorized and GPU-accelerated RMSD and ilRMSD matrix computations.

Provides batched, vectorized Kabsch superposition algorithms using PyTorch
(with CUDA / MPS acceleration) and NumPy fallbacks, avoiding intermediate
file I/O and achieving high speedups on multi-core CPUs and GPUs.
"""

from pathlib import Path

import numpy as np

from haddock import log


def _resolve_torch_device(device: str = "auto") -> str | None:
    """Resolve compute device for PyTorch."""
    try:
        import torch

        if device == "auto":
            if torch.cuda.is_available():
                return "cuda"
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
            return "cpu"
        elif device == "cuda" and torch.cuda.is_available():
            return "cuda"
        elif (
            device == "mps"
            and hasattr(torch.backends, "mps")
            and torch.backends.mps.is_available()
        ):
            return "mps"
        return "cpu"
    except ImportError:
        return None


def compute_rmsd_matrix_torch(
    coords: np.ndarray,
    device: str = "auto",
    batch_size: int = 50000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute pairwise RMSD matrix for N models using batched PyTorch.

    Args:
        coords: NumPy array of shape (N, M, 3) where N is number of models
            and M is number of atoms.
        device: Acceleration device ('auto', 'cuda', 'mps', 'cpu').
        batch_size: Number of structure pairs to process per chunk.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]:
            - i_indices: 1D array of 0-based first structure indices.
            - j_indices: 1D array of 0-based second structure indices.
            - rmsds: 1D array of RMSD values in Angstroms.
    """
    import torch

    target_dev = _resolve_torch_device(device) or "cpu"
    log.info(f"Computing RMSD matrix on device: {target_dev}")

    n_models, n_atoms, _ = coords.shape
    # MPS does not support float64 tensors; fall back to float32 on MPS
    dtype = torch.float32 if str(target_dev).startswith("mps") else torch.float64
    coords_t = torch.as_tensor(coords, dtype=dtype, device=target_dev)

    # Center all structures at origin
    centroids = coords_t.mean(dim=1, keepdim=True)
    coords_c = coords_t - centroids
    # Squared norms for analytic Kabsch formula
    sq_norms = (coords_c**2).sum(dim=(1, 2))  # shape: (N,)

    # Pairwise indices (upper triangular, i < j)
    i_idx, j_idx = torch.triu_indices(n_models, n_models, offset=1, device=target_dev)
    tot_pairs = i_idx.shape[0]

    rmsd_list = []
    for start in range(0, tot_pairs, batch_size):
        end = min(start + batch_size, tot_pairs)
        batch_i = i_idx[start:end]
        batch_j = j_idx[start:end]

        p_i = coords_c[batch_i]  # (B, M, 3)
        p_j = coords_c[batch_j]  # (B, M, 3)

        # Correlation matrices C = p_i^T @ p_j, shape (B, 3, 3)
        c_mat = torch.bmm(p_i.transpose(1, 2), p_j)

        # Batched SVD on (B, 3, 3)
        u, s, vh = torch.linalg.svd(c_mat)
        # Check reflection
        det = torch.det(torch.bmm(u, vh))
        d = torch.sign(det)

        trace = s[:, 0] + s[:, 1] + d * s[:, 2]
        e0 = sq_norms[batch_i] + sq_norms[batch_j]
        diff = torch.clamp((e0 - 2.0 * trace) / float(n_atoms), min=0.0)
        batch_rmsd = torch.sqrt(diff)
        rmsd_list.append(batch_rmsd.cpu())

    all_rmsds = torch.cat(rmsd_list).numpy().astype(np.float64, copy=False)
    return i_idx.cpu().numpy(), j_idx.cpu().numpy(), all_rmsds


def compute_rmsd_matrix_numpy(
    coords: np.ndarray,
    batch_size: int = 10000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute pairwise RMSD matrix for N models using NumPy fallback.

    Args:
        coords: NumPy array of shape (N, M, 3).
        batch_size: Batch size for processing pairs.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]:
            - i_indices: 1D array of 0-based first structure indices.
            - j_indices: 1D array of 0-based second structure indices.
            - rmsds: 1D array of RMSD values in Angstroms.
    """
    n_models, n_atoms, _ = coords.shape
    coords_f64 = coords.astype(np.float64)
    coords_c = coords_f64 - coords_f64.mean(axis=1, keepdims=True)
    sq_norms = (coords_c**2).sum(axis=(1, 2))

    i_idx, j_idx = np.triu_indices(n_models, k=1)
    tot_pairs = len(i_idx)

    rmsds = np.empty(tot_pairs, dtype=np.float64)
    for start in range(0, tot_pairs, batch_size):
        end = min(start + batch_size, tot_pairs)
        b_i = i_idx[start:end]
        b_j = j_idx[start:end]

        p_i = coords_c[b_i]
        p_j = coords_c[b_j]
        # (B, 3, 3)
        c_mat = np.matmul(np.transpose(p_i, (0, 2, 1)), p_j)

        u, s, vh = np.linalg.svd(c_mat)
        det = np.linalg.det(np.matmul(u, vh))
        d = np.sign(det)

        trace = s[:, 0] + s[:, 1] + d * s[:, 2]
        e0 = sq_norms[b_i] + sq_norms[b_j]
        diff = np.maximum(0.0, (e0 - 2.0 * trace) / float(n_atoms))
        rmsds[start:end] = np.sqrt(diff)

    return i_idx, j_idx, rmsds


def compute_rmsd_matrix(
    coords: np.ndarray,
    use_gpu: bool = True,
    device: str = "auto",
    batch_size: int = 50000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate pairwise RMSD matrix for all model pairs.

    Automatically uses GPU tensor operations via PyTorch if available,
    falling back to NumPy.

    Args:
        coords: NumPy array of shape (N, M, 3).
        use_gpu: Whether to utilize GPU acceleration.
        device: Compute device ('auto', 'cuda', 'mps', 'cpu').
        batch_size: Number of pairs per batch chunk.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: (i_indices, j_indices, rmsds)
    """
    if use_gpu:
        try:
            return compute_rmsd_matrix_torch(
                coords, device=device, batch_size=batch_size
            )
        except (ImportError, RuntimeError) as err:
            log.warning(
                f"PyTorch GPU RMSD calculation failed ({err}). Falling back to NumPy."
            )

    return compute_rmsd_matrix_numpy(coords, batch_size=batch_size)


def compute_ilrmsd_matrix_torch(
    rec_coords: np.ndarray,
    lig_coords: np.ndarray,
    device: str = "auto",
    batch_size: int = 25000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute pairwise interface-ligand RMSD (ilRMSD) using batched PyTorch.

    Receptor coordinates are superposed using Kabsch alignment, and the
    resulting optimal rotation and translation are applied to the ligand
    coordinates to measure ligand RMSD.

    Args:
        rec_coords: Receptor coordinates of shape (N, M_rec, 3).
        lig_coords: Ligand coordinates of shape (N, M_lig, 3).
        device: Acceleration device ('auto', 'cuda', 'mps', 'cpu').
        batch_size: Number of pairs to process per chunk.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: (i_indices, j_indices, rmsds)
    """
    import torch

    target_dev = _resolve_torch_device(device) or "cpu"
    log.info(f"Computing ilRMSD matrix on device: {target_dev}")

    n_models, _, _ = rec_coords.shape
    _, n_lig_atoms, _ = lig_coords.shape

    rec_t = torch.as_tensor(rec_coords, dtype=torch.float32, device=target_dev)
    lig_t = torch.as_tensor(lig_coords, dtype=torch.float32, device=target_dev)

    rec_centers = rec_t.mean(dim=1, keepdim=True)
    rec_c = rec_t - rec_centers

    i_idx, j_idx = torch.triu_indices(n_models, n_models, offset=1, device=target_dev)
    tot_pairs = i_idx.shape[0]

    rmsd_list = []
    for start in range(0, tot_pairs, batch_size):
        end = min(start + batch_size, tot_pairs)
        b_i = i_idx[start:end]
        b_j = j_idx[start:end]

        r_i = rec_c[b_i]  # (B, M_rec, 3)
        r_j = rec_c[b_j]  # (B, M_rec, 3)

        # Correlation matrix to align receptor j onto receptor i
        # C = r_j^T @ r_i
        c_mat = torch.bmm(r_j.transpose(1, 2), r_i)
        u, _, vh = torch.linalg.svd(c_mat)
        det = torch.det(torch.bmm(vh.transpose(1, 2), u.transpose(1, 2)))
        d = torch.sign(det)

        # Optimal rotation matrix Rot: (B, 3, 3)
        diag = torch.ones((len(b_i), 3), dtype=torch.float32, device=target_dev)
        diag[:, 2] = d
        diag_mat = torch.diag_embed(diag)
        rot = torch.bmm(torch.bmm(vh.transpose(1, 2), diag_mat), u.transpose(1, 2))

        # Apply rotation and translation to ligand j:
        # l_j_aligned = (l_j - rec_center_j) @ rot^T + rec_center_i
        l_i = lig_t[b_i]
        l_j = lig_t[b_j]
        l_j_c = l_j - rec_centers[b_j]
        l_j_rot = torch.bmm(l_j_c, rot.transpose(1, 2))
        l_j_aligned = l_j_rot + rec_centers[b_i]

        diff_sq = ((l_i - l_j_aligned) ** 2).sum(dim=(1, 2))
        b_rmsd = torch.sqrt(torch.clamp(diff_sq / float(n_lig_atoms), min=0.0))
        rmsd_list.append(b_rmsd.cpu())

    all_rmsds = torch.cat(rmsd_list).numpy()
    return i_idx.cpu().numpy(), j_idx.cpu().numpy(), all_rmsds


def compute_ilrmsd_matrix(
    rec_coords: np.ndarray,
    lig_coords: np.ndarray,
    use_gpu: bool = True,
    device: str = "auto",
    batch_size: int = 25000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate pairwise interface-ligand RMSD matrix for all model pairs.

    Args:
        rec_coords: Receptor coordinates of shape (N, M_rec, 3).
        lig_coords: Ligand coordinates of shape (N, M_lig, 3).
        use_gpu: Whether to utilize GPU acceleration.
        device: Acceleration device ('auto', 'cuda', 'mps', 'cpu').
        batch_size: Number of pairs per batch chunk.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: (i_indices, j_indices, rmsds)
    """
    if use_gpu:
        try:
            return compute_ilrmsd_matrix_torch(
                rec_coords, lig_coords, device=device, batch_size=batch_size
            )
        except (ImportError, RuntimeError) as err:
            log.warning(f"PyTorch GPU ilRMSD calculation failed ({err}).")

    raise RuntimeError("PyTorch is required for GPU-accelerated ilRMSD calculation.")


def write_rmsd_matrix_file(
    output_path: Path | str,
    i_indices: np.ndarray,
    j_indices: np.ndarray,
    rmsds: np.ndarray,
) -> None:
    """Write RMSD matrix in HADDOCK3 standard format (1-based indices).

    Format:
        model1 model2 rmsd_value

    Args:
        output_path: Target output matrix filepath.
        i_indices: 0-based first structure indices.
        j_indices: 0-based second structure indices.
        rmsds: RMSD values in Angstroms.
    """
    with open(output_path, "w") as fh:
        fh.writelines(
            f"{i + 1} {j + 1} {r:.3f}\n" for i, j, r in zip(i_indices, j_indices, rmsds)
        )
