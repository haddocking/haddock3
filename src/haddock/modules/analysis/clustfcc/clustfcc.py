import os

import numpy as np

from haddock import log
from fcc.scripts import cluster_fcc


from pathlib import Path

def iterate_clustering(pool, threshold_param):
    """
    Iterate over the clustering process until a cluster is found.

    Parameters
    ----------
    pool : fcc.Pool
        The pool object containing the fcc matrix.
    
    threshold_param : int
        The threshold parameter to start the clustering process.
    
    Returns
    -------
    clusters : list
        A list of clusters.
    
    threshold : int
        The threshold used to obtain the clusters.
    """
    cluster_check = False
    while not cluster_check:
        for threshold in range(threshold_param, 0, -1):
            log.info(f'Clustering with threshold={threshold}')
            _, clusters = cluster_fcc.cluster_elements(
                pool,
                threshold=threshold,
                )
            if not clusters:
                log.info(
                    "[WARNING] No cluster was found, decreasing threshold!"
                    )
            else:
                cluster_check = True
                # pass the actual threshold back to the param dict
                #  because it will be use in the detailed output
                break
        if not cluster_check:
            # No cluster was obtained in any threshold
            cluster_check = True
    return clusters, threshold

def write_clusters(clusters):
    """
    Write the clusters to the cluster.out file.

    Parameters
    ----------
    clusters : list
        A list of clusters.
    
    Returns
    -------
    None
    """
    # write the classic output file for compatibility reasons
    log.info('Saving output to cluster.out')
    cluster_out = Path('cluster.out')
    with open(cluster_out, 'w') as fh:
        cluster_fcc.output_clusters(fh, clusters)
    fh.close()

def get_cluster_centers(clusters, models_to_cluster):
    """
    Get the cluster centers and the cluster dictionary.

    Parameters
    ----------
    clusters : list
        A list of clusters.
    
    models_to_cluster : list
        A list of models to cluster.
    
    Returns
    -------
    clt_dic : dict
        A dictionary containing the clusters.
    
    clt_centers : dict
        A dictionary containing the cluster centers.
    """
    clt_dic = {}
    clt_centers = {}
    # iterate over the clusters
    for clt in clusters:
        cluster_id = clt.name
        cluster_center_id = clt.center.name - 1
        cluster_center_pdb = models_to_cluster[cluster_center_id]

        clt_dic[cluster_id] = []
        clt_centers[cluster_id] = cluster_center_pdb
        clt_dic[cluster_id].append(cluster_center_pdb)
        # iterate over the models in the cluster
        for model in clt.members:
            model_id = model.name
            model_pdb = models_to_cluster[model_id - 1]
            clt_dic[cluster_id].append(model_pdb)
    return clt_dic, clt_centers


def write_clustfcc_file(clusters, clt_centers, clt_dic, params, sorted_score_dic):
    """
    Write the clustfcc.txt file.

    Parameters
    ----------
    clusters : list
        A list of clusters.
    
    clt_centers : dict
        A dictionary containing the cluster centers.
    
    clt_dic : dict
        A dictionary containing the clusters.
    
    params : dict
        A dictionary containing the clustering parameters.
    
    sorted_score_dic : list
        A list of sorted scores.
    
    Returns
    -------
    None
    """
    # Prepare clustfcc.txt
    output_fname = Path('clustfcc.txt')
    output_str = f'### clustfcc output ###{os.linesep}'
    output_str += os.linesep
    output_str += f'Clustering parameters {os.linesep}'
    output_str += (
        "> contact_distance_cutoff="
        f"{params['contact_distance_cutoff']}A"
        f"{os.linesep}")
    output_str += (
        f"> fraction_cutoff={params['fraction_cutoff']}"
        f"{os.linesep}")
    output_str += f"> threshold={params['threshold']}{os.linesep}"
    output_str += (
        f"> strictness={params['strictness']}{os.linesep}")
    output_str += os.linesep
    output_str += (
        "Note: Models marked with * represent the center of the cluster"
        f"{os.linesep}")
    output_str += (
        f"-----------------------------------------------{os.linesep}")
    output_str += os.linesep
    output_str += f'Total # of clusters: {len(clusters)}{os.linesep}'

    for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
        cluster_id, _ = _e
        center_pdb = clt_centers[cluster_id]
        model_score_l = [(e.score, e) for e in clt_dic[cluster_id]]
        model_score_l.sort()
        # subset_score_l = [e[0] for e in model_score_l][:threshold]
        subset_score_l = [e[0] for e in model_score_l][:params['threshold']]
        top_mean_score = np.mean(subset_score_l)
        top_std = np.std(subset_score_l)
        output_str += (
            f"{os.linesep}"
            "-----------------------------------------------"
            f"{os.linesep}"
            f"Cluster {cluster_rank} (#{cluster_id}, "
            f"n={len(model_score_l)}, "
            f"top{params['threshold']}_avg_score = {top_mean_score:.2f} "
            f"+-{top_std:.2f})"
            f"{os.linesep}")
        output_str += os.linesep
        output_str += f'clt_rank\tmodel_name\tscore{os.linesep}'
        for model_ranking, element in enumerate(model_score_l, start=1):
            score, pdb = element
            if pdb.file_name == center_pdb.file_name:
                output_str += (
                    f"{model_ranking}\t{pdb.file_name}\t{score:.2f}\t*"
                    f"{os.linesep}")
            else:
                output_str += (
                    f"{model_ranking}\t{pdb.file_name}\t{score:.2f}"
                    f"{os.linesep}")
    output_str += (
        "-----------------------------------------------"
        f"{os.linesep}")

    log.info('Saving detailed output to clustfcc.txt')
    with open(output_fname, 'w') as out_fh:
        out_fh.write(output_str)

    return