"""haddock3-re clustrmsd subcommand."""
from pathlib import Path

import numpy as np

from haddock import log
from haddock.core.typing import Union
from haddock.gear.config import load as read_config
from haddock.gear.config import save as save_config
from haddock.libs.libclust import (
    add_cluster_info,
    clustrmsd_tolerance_params,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libinteractive import look_for_capri, rewrite_capri_tables
from haddock.libs.libontology import ModuleIO
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_clusters,
    iterate_min_population,
    write_clusters,
    write_clustrmsd_file,
    )


def add_clustrmsd_arguments(clustrmsd_subcommand):
    """Add arguments to the clustrmsd subcommand."""
    clustrmsd_subcommand.add_argument(
        "clustrmsd_dir",
        help="The clustrmsd directory to recluster.",
        )
    
    clustrmsd_subcommand.add_argument(
        "-n",
        "--n_clusters",
        help="number of clusters to generate.",
        required=False,
        type=int,
        )
    
    clustrmsd_subcommand.add_argument(
        "-d",
        "--clust_cutoff",
        help="clustering cutoff distance.",
        required=False,
        type=float,
        )
    
    clustrmsd_subcommand.add_argument(
        "-t",
        "--min_population",
        help="minimum cluster population.",
        required=False,
        type=int,
        )

    return clustrmsd_subcommand


def reclustrmsd(
        clustrmsd_dir: str,
        n_clusters: Union[bool, int] = None,
        clust_cutoff: Union[bool, float] = None,
        min_population: Union[bool, int] = None,
        ) -> Path:
    """
    Recluster the models in the clustrmsd directory.
    
    Parameters
    ----------
    clustrmsd_dir : str
        Path to the clustrmsd directory.
    
    n_clusters : Union[bool, int]
        Number of clusters to generate.
    
    clust_cutoff : Union[bool, float]
        Clustering cutoff distance.
    
    min_population : Union[bool, int]
        Cluster population min_population.
    
    Returns
    -------
    outdir : Path
        Path to the interactive directory.
    """
    log.info(f"Reclustering {clustrmsd_dir}")

    run_dir = Path(clustrmsd_dir).parent
    clustrmsd_name = Path(clustrmsd_dir).name
    # create the interactive folder
    outdir = Path(run_dir, f"{clustrmsd_name}_interactive")
    outdir.mkdir(exist_ok=True)

    # create an io object
    io = ModuleIO()
    filename = Path(clustrmsd_dir, "io.json")
    io.load(filename)
    models = io.input

    # load the original clustering parameters via json
    clustrmsd_params = read_config(Path(clustrmsd_dir, "params.cfg"))
    key = list(clustrmsd_params['final_cfg'].keys())[0]
    clustrmsd_params = clustrmsd_params['final_cfg'][key]
    log.info(f"Previous clustering parameters: {clustrmsd_params}")

    # setting previous tolerance, just in case no new parameters are given
    tolerance_param_name, tolerance = clustrmsd_tolerance_params(
        clustrmsd_params,
        )

    # adjust the parameters
    if n_clusters is not None:
        clustrmsd_params["n_clusters"] = n_clusters
        clustrmsd_params["criterion"] = "maxclust"
        tolerance = n_clusters
    else:
        if clust_cutoff is not None:
            clustrmsd_params["clust_cutoff"] = clust_cutoff
            clustrmsd_params["criterion"] = "distance"
            tolerance = clust_cutoff
    
    if min_population is not None:
        clustrmsd_params["min_population"] = min_population

    log.info(
        f"Clustering with {tolerance_param_name} = {tolerance}, "
        f"and criterion {clustrmsd_params['criterion']}"
        )
    
    # load the clustering dendrogram
    dendrogram = np.loadtxt(Path(clustrmsd_dir, "dendrogram.txt"))

    # get the clusters
    cluster_arr = get_clusters(
        dendrogram,
        tolerance,
        clustrmsd_params["criterion"],
        )
    log.info(f"clusters {cluster_arr}")

    if clustrmsd_params['criterion'] == "distance":
        cluster_arr, min_population = iterate_min_population(
            cluster_arr,
            clustrmsd_params['min_population']
            )
        clustrmsd_params['min_population'] = min_population
    log.info(f"Updated clustering parameters = {clustrmsd_params}")
    
    # processing the clusters
    unq_clusters = np.unique(cluster_arr)  # contains -1 (unclustered)
    clusters = [c for c in unq_clusters if c != -1]
    log.info(f"clusters = {clusters}")

    clt_dic, cluster_centers = write_clusters(
        clusters,
        cluster_arr,
        models,
        out_filename=Path(outdir, "cluster.out"),
        rmsd_matrix=None,
        centers=False
        )

    score_dic, sorted_score_dic = rank_clusters(
        clt_dic,
        clustrmsd_params["min_population"]
        )
    
    output_models = add_cluster_info(sorted_score_dic, clt_dic)
    
    write_structure_list(
        models,
        output_models,
        out_fname=Path(outdir, "clustrmsd.tsv"),
        )
    
    write_clustrmsd_file(
        clusters,
        clt_dic,
        cluster_centers,
        score_dic,
        sorted_score_dic,
        clustrmsd_params,
        output_fname=Path(outdir, "clustrmsd.txt"),
        )

    # save the io.json file
    io.save(outdir)

    # save the updated parameters in a json file
    save_config(clustrmsd_params, Path(outdir, "params.cfg"))

    # analysis
    clustrmsd_id = int(clustrmsd_name.split("_")[0])
    caprieval_folder = look_for_capri(run_dir, clustrmsd_id)
    if caprieval_folder:
        log.info("Rewriting capri tables")
        rewrite_capri_tables(caprieval_folder, clt_dic, outdir)

    return outdir
