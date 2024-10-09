"""
Library of functions related to the clustering modules.

Main functions
--------------

* :py:func:`write_unclustered_list`
* :py:func:`plot_cluster_matrix`
"""

import os
from pathlib import Path

from haddock import log
from haddock.core.typing import FilePath, Union, ParamDictT, Optional
from haddock.libs.libontology import PDBFile
from haddock.libs.libplots import heatmap_plotly

import numpy as np
from scipy.spatial.distance import squareform


MAX_NB_ENTRY_HTML_MATRIX = 3100


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


def plot_cluster_matrix(
        matrix_path: Union[Path, FilePath, str],
        final_order_idx: list[int],
        labels: list[str],
        dttype: str = '',
        diag_fill: Union[int, float] = 1,
        color_scale: str = "Blues",
        reverse: bool = False,
        output_fname: Union[str, Path, FilePath] = 'clust_matrix',
        matrix_cluster_dt: Optional[list[list[list[int]]]] = None,
        cluster_limits: Optional[list[dict[str, float]]] = None,
        ) -> Optional[str]:
    """Plot a plotly heatmap of a matrix file.

    Parameters
    ----------
    matrix_path : Union[Path, FilePath, str]
        Path to a half-matrix
    final_order_idx : list[int]
        Index orders
    labels : list[str]
        Ordered labels
    dttype : str
        Name of the data type, by default ``
    color_scale : str, optional
        Color scale for the plot, by default "Blues"
    reversed : bool, optional
        Should the color scale be reversed ?, by default False
    output_fname : Union[str, Path, FilePath], optional
        Name of the output file to generate, by default 'clust_matrix.html'
    matrix_cluster_dt: Optional[list[list[list[int]]]]
        A matrix of cluster ids, used for extra hover annotation in plotly.
    cluster_limits: Optional[list[dict[str, float]]]
        A list of dict enabling to draw lines separating cluster ids.

    Return
    ------
    output_fname_ext : str
        Path to the generated file containing the figure.
    """
    # Check that we will be able to generate a functional interactive plot
    if len(final_order_idx) > MAX_NB_ENTRY_HTML_MATRIX:
        return None

    upper_diag, lower_diag = [], []
    # Read matrix
    with open(matrix_path, 'r') as f:
        # Loop over lines
        for _ in f:
            # Split line
            s_ = _.strip().split()
            # Point first value
            uv = float(s_[2])
            # Point second value (if exists)
            lv = float(s_[3]) if len(s_) == 4 else uv
            # Hold them
            upper_diag.append(uv)
            lower_diag.append(lv)

    # Genereate full matrix from N*(N-1)/2 vector
    upper_matrix = squareform(upper_diag)
    lower_matrix = squareform(lower_diag)
    # Update diagonal with data
    np.fill_diagonal(upper_matrix, diag_fill)

    # Full matrix (lower triangle + upper triangle)
    full_matrix = np.tril(lower_matrix, k=-1) + np.triu(upper_matrix)

    # Extract submatrix of selected models and re-order them
    submat = full_matrix[np.ix_(final_order_idx, final_order_idx)]

    # Check if must reverse the colorscale
    if reverse:
        if color_scale[-2:] == '_r':
            color_scale = color_scale[:-2]
        else:
            color_scale += '_r'

    # Define hovering tempalte string
    if matrix_cluster_dt:
        hovertemplate = (
            f'                                {dttype}: %{{z}} <br>'
            f' Model1: %{{x}}    ClusterID: %{{customdata[0]}} <br>'
            f' Model2: %{{y}}    ClusterID: %{{customdata[1]}} '
            '<extra></extra>'
            )
    else:
        hovertemplate = (
            f'                       {dttype}: %{{z}} <br>'
            f' Model1: %{{x}} <br>'
            f' Model2: %{{y}} '
            '<extra></extra>'
            )

    # Generate file name
    output_fname_ext = f"{output_fname}.html"
    # Draw heatmap
    heatmap_plotly(
        submat,
        labels={'color': dttype},
        xlabels=labels,
        ylabels=labels,
        color_scale=color_scale,
        title=f"{dttype} clustering matrix",
        output_fname=output_fname_ext,
        hovertemplate=hovertemplate,
        customdata=matrix_cluster_dt,
        delineation_traces=cluster_limits,
        )
    # Return generated filepath
    return output_fname_ext


def get_cluster_matrix_plot_clt_dt(
        cluster_ids: list[int],
        ) -> tuple[list[list[list[int]]], list[dict[str, float]]]:
    """Generate cluster matrix data for plotly.

    Parameters
    ----------
    cluster_ids : list[int]
        List containing ordered cluster ids.

    Returns
    -------
    matrix_cluster_dt: list[list[list[int]]]
        A matrix of cluster ids, used for plotly.
    
    cluster_limits: list[dict[str, float]]]
        Boundaries to draw lines between clusters with plotly.
    """
    # Set custom data
    matrix_cluster_dt = [
        [[clix, cliy] for clix in cluster_ids]
        for cliy in cluster_ids
        ]
    # Build delineation lines
    del_ind = -0.5
    del_posi = []
    current_clid = cluster_ids[0]
    for clid in cluster_ids:
        if clid != current_clid:
            del_posi.append(del_ind)
            current_clid = clid
        del_ind += 1
    cluster_limits = [
        {
            "x0": delpos,
            "x1": delpos,
            "y0": -0.5,
            "y1": len(cluster_ids) - 0.5,
            }
        for delpos in del_posi
        ] + [
            {
                "y0": delpos,
                "y1": delpos,
                "x0": -0.5,
                "x1": len(cluster_ids) - 0.5,
                }
            for delpos in del_posi
        ]
    return matrix_cluster_dt, cluster_limits


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
        parameters: ParamDictT,
        ) -> tuple[str, Union[int, float]]:
    """Provide parameters of interest for clust rmsd.

    Parameters
    ----------
    parameters : ParamDictT
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
