"""RMSD clustering."""
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

from haddock import log
from haddock.libs.libontology import RMSDFile


def read_matrix(rmsd_matrix):
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
    with open(filename, 'r') as mf:
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
    log.info('Clustering dendrogram...')
    cluster_list = fcluster(dendrogram, t=tolerance, criterion=criterion)
    return cluster_list


def cond_index(i, j, n):
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


def get_cluster_center(npw, n_obs, rmsd_matrix):
    """
    Get the cluster centers.

    Parameters
    ----------
    clusters : list
        List of clusters.
    cluster_list : list
        List of clustered observations.
    rmsd_matrix : array
        List of RMSD values.

    Returns
    -------
    cluster_center : int
        Index of cluster center
    """
    intra_cl_distances = {el: 0.0 for el in npw}
    
    # iterating over the elements of the cluster
    for m_idx in range(len(npw)):
        npws = npw[m_idx + 1:]
        pairs = [int(cond_index(npw[m_idx], npw_el, n_obs)) for npw_el in npws]
        for pair_idx in range(len(pairs)):
            intra_cl_distances[npw[m_idx]] += rmsd_matrix[pairs[pair_idx]]
            intra_cl_distances[npws[pair_idx]] += rmsd_matrix[pairs[pair_idx]]
    cluster_center = min(intra_cl_distances, key=intra_cl_distances.get)
    return cluster_center
