"""RMSD clustering."""
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

from haddock import log
from haddock.libs.libontology import RMSDFile


def read_matrix(rmsd_matrix: RMSDFile) -> np.ndarray:
    """
    Read the RMSD matrix.

    Parameters
    ----------
    rmsd_matrix : :obj:`RMSDFile`
        RMSDFile object with the path to the RMSD matrix.

    Returns
    -------
    matrix : :obj:`numpy.ndarray`
        Numpy array with the RMSD matrix.
    """
    if not isinstance(rmsd_matrix, RMSDFile):
        err = f"{type(rmsd_matrix)} is not a RMSDFile object."
        raise TypeError(err)
    filename = Path(rmsd_matrix.path, rmsd_matrix.file_name)
    # count lines
    nlines = sum(1 for line in open(filename))
    log.info(f"input rmsd matrix has {nlines} entries")
    # must be a 1D condensed distance matrix
    d = int(np.ceil(np.sqrt(nlines * 2)))
    if (d * (d - 1) / 2) != nlines:
        err = f"{nlines} is not a valid binomial coefficient"
        raise ValueError(err)
    if nlines != rmsd_matrix.npairs:
        err = f"number of pairs {nlines} != expected ({rmsd_matrix.npairs})"
        raise ValueError(err)
    # creating and filling matrix obj
    matrix = np.zeros((nlines))
    c = 0
    with open(filename, "r") as mf:
        for line in mf:
            data = line.split()
            if len(data) != 3:
                raise ValueError(f"line {line} malformed")
            else:
                matrix[[c]] = float(data[2])
                c += 1
    return matrix


def get_dendrogram(rmsd_matrix, linkage_type):
    """Get the dendrogram."""
    Z = linkage(rmsd_matrix, linkage_type)
    return Z


def get_clusters(dendrogram, tolerance, criterion):
    """Obtain the clusters."""
    log.info("Clustering dendrogram...")
    cluster_arr = fcluster(dendrogram, t=tolerance, criterion=criterion)
    return cluster_arr


def apply_threshold(cluster_arr: np.ndarray, threshold: int) -> np.ndarray:
    """
    Apply threshold to cluster list.

    Parameters
    ----------
    cluster_arr : np.ndarray
        Array of clusters.
    threshold : int
        Threshold value on cluster population.

    Returns
    -------
    cluster_arr : np.ndarray
        Array of clusters (unclustered structures are labelled with -1)
    """
    new_cluster_arr = cluster_arr.copy()
    log.info(f"Applying threshold {threshold} to cluster list")
    cluster_pops = np.unique(cluster_arr, return_counts=True)
    uncl_idx = np.where(cluster_pops[1] < threshold)[0]
    invalid_clusters = cluster_pops[0][uncl_idx]
    log.info(f"Invalid clusters: {invalid_clusters}")
    # replacing invalid clusters with -1
    uncl_models = 0
    for cl_idx, cl_id in enumerate(new_cluster_arr):
        if cl_id in invalid_clusters:
            new_cluster_arr[cl_idx] = -1
            uncl_models += 1
    log.info(f"Threshold applied, {uncl_models} models left unclustered")
    return new_cluster_arr


def iterate_threshold(cluster_arr: np.ndarray, threshold: int) -> np.ndarray:
    """
    Iterate over the threshold values until we find at least one valid cluster.

    Parameters
    ----------
    cluster_arr : np.ndarray
        Array of clusters.
    threshold : int
        Threshold value on cluster population.

    Returns
    -------
    new_cluster_arr : np.ndarray
        Array of clusters (unclustered structures are labelled with -1)
    """
    new_cluster_arr: np.ndarray = np.ndarray([])
    for curr_thr in range(threshold, 0, -1):
        log.info(f"Clustering with threshold={curr_thr}")
        new_cluster_arr = apply_threshold(cluster_arr, curr_thr)
        ncl = len(np.unique(new_cluster_arr))
        if ncl <= 1:  # contains -1 (unclustered)
            log.warning(f"No clusters found with threshold={curr_thr}")
        else:
            break
    return new_cluster_arr


def cond_index(i: int, j: int, n: int) -> float:
    """
    Get the condensed index from two matrix indexes.

    Parameters
    ----------
    i : int
        Index of the first element.
    j : int
        Index of the second element.
    n : int
        Number of observations.
    """
    return n * (n - 1) / 2 - (n - i) * (n - i - 1) / 2 + j - i - 1


def get_cluster_center(npw: np.ndarray, n_obs: int, rmsd_matrix: np.ndarray) -> int:
    """
    Get the cluster centers.

    Parameters
    ----------
    npw: np.ndarray
        Indexes of the cluster over cluster_list array
    n_obs : int
        Number of overall observations (models).
    rmsd_matrix : np.ndarray
        RMSD matrix.

    Returns
    -------
    cluster_center : int
        Index of cluster center
    """
    intra_cl_distances = {el: 0.0 for el in npw}

    # iterating over the elements of the cluster
    for m_idx in range(len(npw)):
        npws = npw[m_idx + 1 :]
        pairs = [int(cond_index(npw[m_idx], npw_el, n_obs)) for npw_el in npws]
        for pair_idx in range(len(pairs)):
            intra_cl_distances[npw[m_idx]] += rmsd_matrix[pairs[pair_idx]]
            intra_cl_distances[npws[pair_idx]] += rmsd_matrix[pairs[pair_idx]]
    cluster_center = min(intra_cl_distances, key=intra_cl_distances.get)  # type: ignore
    return cluster_center
