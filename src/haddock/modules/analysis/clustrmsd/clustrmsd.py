"""RMSD clustering."""
import os
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

from haddock import log
from haddock.libs.libontology import RMSDFile


def get_matrix_path(rmsd_matrix: RMSDFile) -> Path:
    """From an RMSDFile object returns the rmsd matrix path.

    Parameters
    ----------
    rmsd_matrix : :obj:`RMSDFile`
        RMSDFile object with the path to the RMSD matrix.

    Returns
    -------
    matrix_fpath : Path
        Path to the RMSD matrix

    Raises
    ------
    TypeError
        If input rmsd_matrix is not of type RMSDFile
    """
    if not isinstance(rmsd_matrix, RMSDFile):
        err = f"{type(rmsd_matrix)} is not a RMSDFile object."
        raise TypeError(err)
    matrix_fpath = Path(rmsd_matrix.path, rmsd_matrix.file_name)
    return matrix_fpath


def write_clusters(clusters, cluster_arr, models, rmsd_matrix, out_filename="cluster.out", centers=False):  # noqa: E501
    """
    Write the clusters to a file.

    Parameters
    ----------
    clusters : list
        List of clusters to write.
    cluster_arr : np.array
        Array with the cluster assignment for each model.
    models : list
        List of models.
    rmsd_matrix : np.array
        RMSD matrix.
    out_filename : str, optional
        Output filename. The default is "cluster.out".
    centers : bool, optional
        Whether to calculate the cluster centers. The default is False.
    
    Returns
    -------
    clt_dic : :obj:`dict`
        Dictionary with the clusters.
    
    cluster_centers : :obj:`dict`
        Dictionary with the cluster ID as key and the cluster center as value.
    """
    n_obs = len(cluster_arr)
    cluster_centers = {}
    # preparing output
    clt_dic = {}
    log.info(f'Saving output to {out_filename}')
    cluster_out = Path(out_filename)
    with open(cluster_out, 'w') as fh:
        for cl_id in clusters:
            if cl_id != -1:
                npw = np.where(cluster_arr == cl_id)[0]
                clt_dic[cl_id] = [models[n] for n in npw]
                fh.write(f"Cluster {cl_id} -> ")

                # find the cluster center
                if centers:
                    clt_center = get_cluster_center(npw, n_obs, rmsd_matrix)
                    cluster_centers[cl_id] = models[clt_center].file_name
                
                # write the cluster
                for el in npw[:-1]:
                    fh.write(f"{el + 1} ")
                fh.write(f"{npw[-1] + 1}")
                fh.write(os.linesep)

    return clt_dic, cluster_centers


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
    filename = get_matrix_path(rmsd_matrix)
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
    return np.array(matrix)


def get_dendrogram(rmsd_matrix, linkage_type):
    """Get and save the dendrogram.
    
    Parameters
    ----------
    rmsd_matrix : :obj:`numpy.ndarray`
        Numpy array with the RMSD matrix.
    
    linkage_type : str
        Linkage type for the clustering.
    
    Returns
    -------
    Z : :obj:`numpy.ndarray`
        Numpy array with the dendrogram.
    """
    Z = linkage(rmsd_matrix, linkage_type)
    np.savetxt("dendrogram.txt", Z, fmt='%.5f')
    return Z


def get_clusters(dendrogram, tolerance, criterion):
    """Obtain the clusters."""
    log.info("Clustering dendrogram...")
    cluster_arr = fcluster(dendrogram, t=tolerance, criterion=criterion)
    return cluster_arr


def apply_min_population(
        cluster_arr: np.ndarray,
        min_population: int,
        ) -> np.ndarray:
    """
    Apply min_population to cluster list.

    Parameters
    ----------
    cluster_arr : np.ndarray
        Array of clusters.
    min_population : int
        min_population value on cluster population.

    Returns
    -------
    cluster_arr : np.ndarray
        Array of clusters (unclustered structures are labelled with -1)
    """
    new_cluster_arr = cluster_arr.copy()
    log.info(f"Applying min_population {min_population} to cluster list")
    cluster_pops = np.unique(cluster_arr, return_counts=True)
    uncl_idx = np.where(cluster_pops[1] < min_population)[0]
    invalid_clusters = cluster_pops[0][uncl_idx]
    log.info(f"Invalid clusters: {invalid_clusters}")
    # replacing invalid clusters with -1
    uncl_models = 0
    for cl_idx, cl_id in enumerate(new_cluster_arr):
        if cl_id in invalid_clusters:
            new_cluster_arr[cl_idx] = -1
            uncl_models += 1
    log.info(f"min_population applied, {uncl_models} models left unclustered")
    return new_cluster_arr


def iterate_min_population(
        cluster_arr: np.ndarray,
        min_population: int,
        ) -> np.ndarray:
    """
    Find one valid valuster satisfying the min_population parameter.

    Logic: Iterate over the min_population values until we find
    at least one valid cluster.

    Parameters
    ----------
    cluster_arr : np.ndarray
        Array of clusters.
    min_population : int
        min_population value on cluster population.

    Returns
    -------
    new_cluster_arr : np.ndarray
        Array of clusters (unclustered structures are labelled with -1)
    """
    new_cluster_arr: np.ndarray = np.ndarray([])
    for curr_thr in range(min_population, 0, -1):
        log.info(f"Clustering with min_population={curr_thr}")
        new_cluster_arr = apply_min_population(cluster_arr, curr_thr)
        ncl = len(np.unique(new_cluster_arr))
        if -1 in new_cluster_arr:  # contains -1 (unclustered)
            ncl -= 1
        # if no clusters found, try with a lower min_population
        if ncl == 0:
            log.warning(f'No clusters found with min_population={curr_thr}')
        else:
            break
    return new_cluster_arr, curr_thr


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


def get_cluster_center(
        npw: np.ndarray,
        n_obs: int,
        rmsd_matrix: np.ndarray,
        ) -> int:
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
        npws = npw[m_idx + 1:]
        pairs = [int(cond_index(npw[m_idx], npw_el, n_obs)) for npw_el in npws]
        for pair_idx in range(len(pairs)):
            intra_cl_distances[npw[m_idx]] += rmsd_matrix[pairs[pair_idx]]
            intra_cl_distances[npws[pair_idx]] += rmsd_matrix[pairs[pair_idx]]
    cluster_center = min(intra_cl_distances, key=intra_cl_distances.get)  # type: ignore  # noqa : E501
    return cluster_center


def write_clustrmsd_file(clusters, clt_dic, cluster_centers, score_dic, sorted_score_dic, params, output_fname='clustrmsd.txt'):  # noqa E501
    """
    Write the clustrmsd.txt file.

    Parameters
    ----------
    clusters : np.ndarray
        Array of clusters.
    
    clt_dic : dict
        Dictionary with the clusters.
    
    cluster_centers : dict
        Dictionary with the cluster centers.
    
    score_dic : dict
        Dictionary with the scores.
    
    sorted_score_dic : dict
        Dictionary with the sorted scores.
    
    params : dict
        Dictionary with the clustering parameters.
    
    output_fname : str
        Output filename.
    """
    # Prepare clustrmsd.txt
    output_str = f'### clustrmsd output ###{os.linesep}'
    output_str += os.linesep
    output_str += f'Clustering parameters {os.linesep}'
    output_str += f"> linkage_type={params['linkage']}{os.linesep}"
    output_str += f"> criterion={params['criterion']}{os.linesep}"
    if params['criterion'] == "distance":
        output_str += f"> clust_cutoff={params['clust_cutoff']:.2f}{os.linesep}"
    else:
        output_str += f"> n_clusters={params['n_clusters']}{os.linesep}"
    output_str += f"> min_population={params['min_population']}{os.linesep}"
    output_str += os.linesep
    
    output_str += (
        f"-----------------------------------------------{os.linesep}")
    output_str += os.linesep
    output_str += f'Total # of clusters: {len(clusters)}{os.linesep}'
    for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
        cluster_id, _ = _e
        
        model_score_l = [(e.score, e) for e in clt_dic[cluster_id]]
        model_score_l.sort()
        top_score = score_dic[cluster_id]

        output_str += (
            f"{os.linesep}"
            "-----------------------------------------------"
            f"{os.linesep}"
            f"Cluster {cluster_rank} (#{cluster_id}, "
            f"n={len(model_score_l)}, "
            f"top{params['min_population']}_avg_score = {top_score:.2f})"
            f"{os.linesep}")
        output_str += os.linesep
        output_str += f'clt_rank\tmodel_name\tscore{os.linesep}'
        for model_ranking, element in enumerate(model_score_l, start=1):
            score, pdb = element
            output_str += (
                f"{model_ranking}\t{pdb.file_name}\t{score:.2f}"
                )
            if cluster_centers:
                # is the model the cluster center?
                if pdb.file_name == cluster_centers[cluster_id]:
                    output_str += "\t*"
            output_str += (f"{os.linesep}")
    output_str += (
        "-----------------------------------------------"
        f"{os.linesep}")
    log.info('Saving detailed output to clustrmsd.txt')
    with open(output_fname, 'w') as out_fh:
        out_fh.write(output_str)


def order_clusters(cluster_arr):
    """
    Order the clusters by population.
    
    The most populated cluster will be assigned the ID 1, the second most
     populated the ID 2, and so on.

    Parameters
    ----------
    cluster_arr : np.ndarray
        Array of clusters.
    
    Returns
    -------
    clusters : list
        List of clusters.
    
    cluster_arr : np.ndarray
        Array of clusters.
    """
    unique_clusters, cluster_counts = np.unique(cluster_arr, return_counts=True)
    
    sorted_indices = np.argsort(-cluster_counts)  # must use negative to sort ascending
    sorted_clusters = unique_clusters[sorted_indices]
    
    # delete -1 from sorted_clusters if present
    sorted_clusters = sorted_clusters[sorted_clusters != -1]
    clusters = []
    # for every element of sorted_clusters I want to assign a new ID
    # to the elements of cluster_arr that match the order of sorted_clusters
    index_dict = {}
    for c in sorted_clusters:
        # index of cluster_arr where the cluster is equal to c
        idx = np.where(cluster_arr == c)
        index_dict[c] = idx
    # now the assignment
    for i, c in enumerate(sorted_clusters):
        clusters.append(i + 1)
        cluster_arr[index_dict[c]] = i + 1
    return clusters, cluster_arr
