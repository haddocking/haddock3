"""RMSD clustering."""
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

from haddock import log
from haddock.libs.libontology import RMSDFile


def read_matrix(rmsd_matrix):
    """Read the RMSD matrix."""
    if not isinstance(rmsd_matrix, RMSDFile):
        err = f"{type(rmsd_matrix)} is not a RMSDFile object."
        raise Exception(err)
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
                raise Exception(f"line {line} malformed")
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
