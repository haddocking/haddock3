"""Set of functions related to the selection of top clusters."""

import os
import math
from pathlib import Path
from haddock.core.typing import Union
from haddock.libs.libontology import PDBFile


def select_top_clusts_models(
        sortby: str,
        models_to_select: list[PDBFile],
        top_clusters: int,
        top_models: Union[int, float],
        ) -> tuple[list[PDBFile], list[str]]:
    """Select best clusters based on structures scores.

    Parameters
    ----------
    sortby : str
        How to order clusters: by `score` or by `size`.
    models_to_select : list[PDBFile]
        List of input models on which selection must be performed.
    top_clusters : int
        Number of best clusters to take into account.
    top_models : int
        Number of best models in each cluster to take into account.

    Returns
    -------
    models_to_export : list[PDBFile]
        List of PDBfiles to export.
    notes : list[str]
        List of notes to be printed.
    """
    notes: list[str] = []
    by_clusters = map_clusters_models(models_to_select)

    # Get cluster order
    if sortby == "size":
        cluster_rankings = size_clust_order(by_clusters)
    else:
        cluster_rankings = rank_clust_order(by_clusters)

    # Check if number of clusters >= set of rank
    if top_clusters >= len(cluster_rankings):
        # select all clusters
        cluster_rankins_str = ",".join(map(str, cluster_rankings))
        notes.append(f"Selecting all clusters: {cluster_rankins_str}")
    else:
        # select top_clusters clusters
        cluster_rankings = cluster_rankings[:top_clusters]
        cluster_rankins_str = ",".join(map(str, cluster_rankings))
        notes.append(
            f"Selecting top {top_clusters} clusters: "
            f"{cluster_rankins_str}"
            )

    # Initiate set of selected models to export
    models_to_export: list[PDBFile] = []
    # Loop over cluster ranks
    for clt_rank in cluster_rankings:
        # Sort models by model rank
        clt_mdls, note = sort_models(by_clusters[clt_rank])
        if note:
            notes.append(note)

        # Set new ranks to models
        for mdl_rank, pdb in enumerate(clt_mdls, start=1):
            pdb.clt_rank = clt_rank
            pdb.clt_model_rank = mdl_rank
        # In case number of models is not set (nan.)
        if math.isnan(top_models):
            for pdb in clt_mdls:
                models_to_export.append(pdb)
        # In case number of models is a integer
        else:
            # Loop over first `top_models` models
            for pdb in clt_mdls[:top_models]:
                models_to_export.append(pdb)
    return models_to_export, notes


def sort_models(
        models: list[PDBFile]
        ) -> tuple[list[PDBFile], Union[None, str]]:
    """Sort models based on their rank in cluster.

    Parameters
    ----------
    models : list[PDBFile]
        List of input models on which ordering must be performed.

    Returns
    -------
    sorted_mdls : list[PDBFile]
        List of sorted models.
    """
    note: Union[None, str] = None
    try:
        sorted_mdls = sorted(
            models,
            key=lambda k: k.clt_model_rank,
            )
    except TypeError:
        note = 'model rank unavailable, falling back to input order'
        sorted_mdls = models
    return sorted_mdls, note
    

def rank_clust_order(
        by_clusters: dict[int, list[PDBFile]],
        ) -> list[int]:
    """Select best clusters based on structures scores.

    Parameters
    ----------
    models_to_select : list[PDBFile]
        List of input models on which selection must be performed.
    top_clusters : int
        Number of best clusters to take into account.
    top_models : int
        Number of best models in each cluster to take into account.

    Returns
    -------
    models_to_export : list[PDBFile]
        List of PDBfiles to export.
    notes : list[str]
        List of notes to be printed.
    """
    # Generate set of all cluster rank available
    cluster_rankings = sorted(by_clusters)
    return cluster_rankings


def size_clust_order(
        by_clusters: dict[int, list[PDBFile]],
        ) -> list[int]:
    """Select best clusters based on structures scores.

    Parameters
    ----------
    models_to_select : list[PDBFile]
        List of input models on which selection must be performed.
    top_clusters : int
        Number of best clusters to take into account.
    top_models : int
        Number of best models in each cluster to take into account.

    Returns
    -------
    models_to_export : list[PDBFile]
        List of PDBfiles to export.
    notes : list[str]
        List of notes to be printed.
    """
    # Generate set of all cluster rank available
    cluster_rankings = sorted(
        by_clusters,
        key=lambda k: len(by_clusters[k]),
        reverse=True,
        )
    return cluster_rankings


def map_clusters_models(models: list[PDBFile]) -> dict[int, list[PDBFile]]:
    """Group models by clusters.

    Parameters
    ----------
    models : list[PDBFile]
        List of PDBfiles models to be grouped.

    Returns
    -------
    by_clusters : dict[int, list[PDBFile]]
        _description_
    """
    # Preset dictionary keys
    by_clusters: dict[int, list[PDBFile]] = {
        clrank: []
        for clrank in list(set([pdb.clt_rank for pdb in models]))
        }
    # Loop over models
    for pdb in models:
        # Add model to cluster
        by_clusters[pdb.clt_rank].append(pdb)
    return by_clusters


def write_selected_models(
        output_path: Union[str, Path],
        models: list[PDBFile],
        module_path: Union[str, Path],
        ) -> list[PDBFile]:
    """Dump selected models and new names in a file.

    Parameters
    ----------
    output_path : Union[str, Path]
        Name of tne file to create.
    models : list[PDBFile]
        List of PDBfiles of selected models.
    module_path : Union[str, Path]
        Path of the module.

    Returns
    -------
    models : list[PDBFile]
        Updated list of selected models.
    """
    # dump the models to disk and change their attributes
    with open(output_path, 'w') as fh:
        fh.write("rel_path\tori_name\tcluster_name\tmd5" + os.linesep)
        for model in models:
            name = (
                f"cluster_{model.clt_rank}_model"
                f"_{model.clt_model_rank}.pdb"
                )
            # writing name
            fh.write(
                f"{model.rel_path}\t"
                f"{model.ori_name}\t"
                f"{name}\t"
                f"{model.md5}" + os.linesep
                )
            # changing attributes
            name_path = Path(name)
            name_path.write_text(model.rel_path.read_text())
            model.ori_name = model.file_name
            model.file_name = name
            model.full_name = name
            model.rel_path = Path('..', Path(module_path).name, name)
            model.path = str(Path(".").resolve())
    return models
