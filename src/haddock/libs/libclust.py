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
from haddock.core.typing import FilePath, Union
from haddock.libs.libontology import PDBFile
from haddock.libs.libplots import heatmap_plotly

import numpy as np
from scipy.spatial.distance import squareform


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
        output_fname: Union[str, Path, FilePath] = 'clust_matrix.html',
        ) -> None:
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
    """
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

    # Extract submatrix of selected models and re-order it
    submat = full_matrix[np.ix_(final_order_idx, final_order_idx)]

    # Check if must reverse the colorscale
    if reverse:
        if color_scale[-2:] == '_r':
            color_scale = color_scale[:-2]
        else:
            color_scale += '_r'

    # Draw heatmap
    heatmap_plotly(
        submat,
        labels={'color': dttype},
        xlabels=labels,
        ylabels=labels,
        color_scale=color_scale,
        title=f"{dttype} clustering matrix",
        output_fname=output_fname,
        )
