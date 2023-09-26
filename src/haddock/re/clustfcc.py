import argparse

import numpy as np
import json

from haddock import log
from pathlib import Path

from haddock.gear.config import load as read_config
from fcc.scripts import calc_fcc_matrix, cluster_fcc

from haddock.libs.libontology import ModuleIO

from haddock.libs.libclust import write_structure_list


def add_clustfcc_arguments(clustfcc_subcommand):
    clustfcc_subcommand.add_argument(
    "clustfcc_dir",
    help="The clustfcc directory to recluster.",
    )

    clustfcc_subcommand.add_argument(
        "-f",
        "--fraction",
        help="fraction of common contacts to not be considered a singleton model.",
        required=False,
        type=float,
        )

    clustfcc_subcommand.add_argument(
        "-s",
        "--strictness",
        help="fraction of common contacts to be considered to be part of the same cluster.",
        required=False,
        type=float,
        )
    
    clustfcc_subcommand.add_argument(
        "-t",
        "--threshold",
        help="cluster population threshold.",
        required=False,
        type=int,
        )
    
    return clustfcc_subcommand

def reclustfcc(clustfcc_dir, fraction=None, strictness=None, threshold=None):
    """Recluster the models in the clustfcc directory."""
    log.info(f"Reclustering {clustfcc_dir}")

    # create the interactive folder
    run_dir = Path(clustfcc_dir).parent
    clustfcc_name = Path(clustfcc_dir).name
    outdir = Path(run_dir, f"{clustfcc_name}_interactive")
    outdir.mkdir(exist_ok=True)

    # create an io object
    io = ModuleIO()
    filename = Path(clustfcc_dir, "io.json")
    io.load(filename)
    models = io.input

    # load the original clustering parameters via json
    clustfcc_params = read_config(Path(clustfcc_dir, "params.cfg"))
    key = list(clustfcc_params['final_cfg'].keys())[0]
    clustfcc_params = clustfcc_params['final_cfg'][key]
    log.info(f"Previous clustering parameters: {clustfcc_params}")

    # adjust the parameters
    if fraction is not None:
        clustfcc_params["fraction"] = fraction
    if strictness is not None:
        clustfcc_params["strictness"] = strictness
    if threshold is not None:
        clustfcc_params["threshold"] = threshold
    
    # load the fcc matrix
    pool = cluster_fcc.read_matrix(
            Path(clustfcc_dir, "fcc.matrix"),
            clustfcc_params['fraction_cutoff'],
            clustfcc_params['strictness'],
            )

    cluster_check = False
    while not cluster_check:
        for threshold in range(clustfcc_params['threshold'], 0, -1):
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
                clustfcc_params['threshold'] = threshold
                break
        if not cluster_check:
            # No cluster was obtained in any threshold
            cluster_check = True

    # Prepare output and read the elements
    clt_dic = {}
    if clusters:
        # write the classic output file for compatibility reasons
        log.info('Saving output to cluster.out')
        cluster_out = Path('cluster.out')
        with open(cluster_out, 'w') as fh:
            cluster_fcc.output_clusters(fh, clusters)
        fh.close()
        clt_centers = {}
        for clt in clusters:
            cluster_id = clt.name
            cluster_center_id = clt.center.name - 1
            cluster_center_pdb = models[cluster_center_id]
            clt_dic[cluster_id] = []
            clt_centers[cluster_id] = cluster_center_pdb
            clt_dic[cluster_id].append(cluster_center_pdb)
            for model in clt.members:
                model_id = model.name
                model_pdb = models[model_id - 1]
                clt_dic[cluster_id].append(model_pdb)
        # Rank the clusters
        #  they are sorted by the topX (threshold) models in each cluster
        score_dic = {}
        for clt_id in clt_dic:
            score_l = [p.score for p in clt_dic[clt_id]]
            score_l.sort()
            denom = float(min(threshold, len(score_l)))
            top4_score = sum(score_l[:threshold]) / denom
            score_dic[clt_id] = top4_score
        sorted_score_dic = sorted(score_dic.items(), key=lambda k: k[1])
        # Add this info to the models
        output_models = []
        for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
            cluster_id, _ = _e
            # sort the models by score
            clt_dic[cluster_id].sort()
            # rank the models
            for model_ranking, pdb in enumerate(clt_dic[cluster_id],
                                                start=1):
                pdb.clt_id = cluster_id
                pdb.clt_rank = cluster_rank
                pdb.clt_model_rank = model_ranking
                output_models.append(pdb)
        # Write unclustered structures
        write_structure_list(models,
                             output_models,
                             out_fname=Path(outdir,"clustfcc.tsv"))
    
    return outdir

