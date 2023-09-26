import argparse
import numpy as np
import json

from haddock import log
from pathlib import Path

from haddock.gear.config import load as read_config

from haddock.libs.libontology import ModuleIO

from haddock.modules.analysis.clustrmsd.clustrmsd import get_clusters, rank_clusters, write_clusters
from haddock.libs.libclust import write_structure_list



def add_clustrmsd_arguments(clustrmsd_subcommand):
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

def reclustrmsd(clustrmsd_dir, n_clusters=None, distance=None, threshold=None):
    """Recluster the models in the clustrmsd directory."""
    log.info(f"Reclustering {clustrmsd_dir}")
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

    ## get the clusters
    cluster_arr = get_clusters(dendrogram, clustrmsd_params["tolerance"], clustrmsd_params["criterion"])
    log.info(f"clusters {cluster_arr}")

    # we got the clusters, now we need write down (part of) the information available in the clustrmsd directory
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
    
    # processing the clusters
    unq_clusters = np.unique(cluster_arr)  # contains -1 (unclustered)
    clusters = [c for c in unq_clusters if c != -1]
    log.info(f"clusters = {clusters}")

    clt_dic, cluster_centers = write_clusters(clusters, cluster_arr, models, out_filename=Path(outdir, "cluster.out"), rmsd_matrix = None, centers=False)

    sorted_score_dic = rank_clusters(clt_dic, clustrmsd_params["threshold"])
    
    # add this to the models
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
    
    write_structure_list(models,
                         output_models,
                         out_fname=Path(outdir,"clustrmsd.tsv"))
    # save the io.json file
    io.save(outdir)

    return outdir
    


