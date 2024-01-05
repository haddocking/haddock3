"""
Library of functions related to the clustering modules.

Main functions
--------------

* :py:func:`write_unclustered_list`
"""
import os
from pathlib import Path

from haddock import log
from haddock.core.typing import FilePath, Union
from haddock.libs.libontology import PDBFile


def write_structure_list(input_models: list[PDBFile],
                         clustered_models: list[PDBFile],
                         out_fname: FilePath) -> None:
    """
    Get the list of unclustered structures.
    
    Parameters
    ----------
    input_models : list
        list of input models
    clustered_models : list
        list of clustered models
    """
    output_fname = Path(out_fname)
    output_str = f'rank\tmodel_name\tscore\tcluster_id{os.linesep}'
    structure_list: list[PDBFile] = []
    # checking which input models have not been clustered
    for model in input_models:
        if model not in clustered_models:
            model.clt_id = "-"
            structure_list.append(model)
    # extending and sorting
    structure_list.extend(clustered_models)
    structure_list.sort(key=lambda model: model.score)
    # adding models to output string
    for mdl_rank, mdl in enumerate(structure_list, start=1):
        output_str += (
            f'{mdl_rank}\t{mdl.file_name}\t{mdl.score:.2f}\t{mdl.clt_id}'
            f'{os.linesep}'
            )
    output_str += os.linesep
    log.info(f'Saving structure list to {out_fname}')
    with open(output_fname, 'w') as out_fh:
        out_fh.write(output_str)


def rank_clusters(clt_dic, threshold):
    """
    Rank the clusters by their average score.

    Parameters
    ----------
    clt_dic : :obj:`dict`
        Dictionary with the clusters.
    
    threshold : int
        Number of models to consider for the average score.
    
    Returns
    -------
    score_dic : :obj:`dict`
        Dictionary with the cluster ID as key and the average score as value.
    
    sorted_score_dic : :obj:`list`
        List of tuples with the cluster ID and the average score, sorted by
        the average score.
    """
    score_dic = {}
    for clt_id in clt_dic:
        score_l = [p.score for p in clt_dic[clt_id]]
        score_l.sort()
        denom = float(min(threshold, len(score_l)))
        top4_score = sum(score_l[:threshold]) / denom
        score_dic[clt_id] = top4_score
    
    sorted_score_dic = sorted(score_dic.items(), key=lambda k: k[1])
    return score_dic, sorted_score_dic


def add_cluster_info(sorted_score_dic, clt_dic):
    """
    Add cluster information to the models.

    Parameters
    ----------
    sorted_score_dic : :obj:`list`
        List of tuples with the cluster ID and the average score, sorted by
        the average score.
    
    clt_dic : :obj:`dict`
        Dictionary with the clusters.
    
    Returns
    -------
    output_models : :obj:`list`
        List of models with the cluster information.
    """
    # Add this info to the models
    output_models = []
    for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
        cluster_id, _ = _e
        # sort the models by score
        clt_dic[cluster_id].sort()
        # rank the models
        for model_ranking, pdb in enumerate(clt_dic[cluster_id],
                                            start=1):
            pdb.clt_id = int(cluster_id)
            pdb.clt_rank = cluster_rank
            pdb.clt_model_rank = model_ranking
            output_models.append(pdb)
    return output_models


def clustrmsd_tolerance_params(
        parameters: dict,
        ) -> tuple[str, Union[int, float]]:
    """Provide parameters of interest for clust rmsd.

    Parameters
    ----------
    parameters : dict
        The clustrmsd module parameters

    Returns
    -------
    tuple[str, Union[int, float]]
        Name of the tolerance parameter and its value.
    """
    # adjust the parameters
    if parameters["criterion"] == "maxclust":
        tolerance_param_name = "n_clusters"
        tolerance = parameters[tolerance_param_name]
    else:  # Expected to be parameters["criterion"] == "distance"
        tolerance_param_name = "clust_cutoff"
        tolerance = parameters[tolerance_param_name]
    return tolerance_param_name, tolerance