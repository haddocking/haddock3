"""haddock3-re clustrmsd subcommand."""
from pathlib import Path

import numpy as np

from haddock import log
from haddock.gear.config import load as read_config
from haddock.gear.config import save as save_config
from haddock.libs.libclust import (
    add_cluster_info,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libinteractive import look_for_capri, rewrite_capri_tables
from haddock.libs.libontology import ModuleIO
from haddock.libs.libplots import read_capri_table
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_clusters,
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
        "--distance",
        help="cutoff distance.",
        required=False,
        type=int,
        )
    
    clustrmsd_subcommand.add_argument(
        "-t",
        "--threshold",
        help="cluster population threshold.",
        required=False,
        type=int,
        )

    return clustrmsd_subcommand


def reclustrmsd(clustrmsd_dir, n_clusters=None, distance=None, threshold=None, caprieval_folder=None):
    """
    Recluster the models in the clustrmsd directory.
    
    Parameters
    ----------
    clustrmsd_dir : str
        Path to the clustrmsd directory.
    
    n_clusters : int
        Number of clusters to generate.
    
    distance : int
        Cutoff distance.
    
    threshold : int
        Cluster population threshold.
    
    caprieval_folder : str
        Path to the caprieval folder for quick analysis of the results.
    
    Returns
    -------
    outdir : str
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
    models = io.retrieve_models()

    # load the original clustering parameters via json
    clustrmsd_params = read_config(Path(clustrmsd_dir, "params.cfg"))
    key = list(clustrmsd_params['final_cfg'].keys())[0]
    clustrmsd_params = clustrmsd_params['final_cfg'][key]
    log.info(f"Previous clustering parameters: {clustrmsd_params}")

    # adjust the parameters
    if n_clusters is not None:
        clustrmsd_params["tolerance"] = n_clusters
        clustrmsd_params["criterion"] = "maxclust"
    else:
        if distance is not None:
            clustrmsd_params["tolerance"] = distance
            clustrmsd_params["criterion"] = "distance"

    if threshold is not None:
        clustrmsd_params["threshold"] = threshold
    
    # load the clustering dendrogram
    dendrogram = np.loadtxt(Path(clustrmsd_dir, "dendrogram.txt"))

    # get the clusters
    cluster_arr = get_clusters(
        dendrogram,
        clustrmsd_params["tolerance"],
        clustrmsd_params["criterion"])
    log.info(f"clusters {cluster_arr}")
    
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
        clustrmsd_params["threshold"]
        )
    
    output_models = add_cluster_info(sorted_score_dic, clt_dic)
    
    write_structure_list(models,
                         output_models,
                         out_fname=Path(outdir, "clustrmsd.tsv")
                         )
    
    write_clustrmsd_file(
        clusters,
        clt_dic,
        cluster_centers,
        score_dic,
        sorted_score_dic,
        clustrmsd_params,
        output_fname=Path(outdir, "clustrmsd.txt")
        )

    # save the io.json file
    io.save(outdir)

    # save the updated parameters in a json file
    save_config(clustrmsd_params, Path(outdir, "params.cfg"))

    # analysis
    clustrmsd_id = int(clustrmsd_dir.split("/")[-1].split("_")[0])
    caprieval_folder = look_for_capri(run_dir, clustrmsd_id)
    if caprieval_folder:
        log.info("Rewriting capri tables")
        rewrite_capri_tables(caprieval_folder, clt_dic, outdir)
    


    return outdir
