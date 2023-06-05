"""Cluster modules with FCC."""
import collections
import os
from pathlib import Path

import numpy as np

from haddock import log
from haddock.libs.libclust import write_structure_list
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.clustfcc.clustfcc import (
    calculate_contacts,
    calculate_fcc_matrix,
    cluster_fcc,
)

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with FCC."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if FCC is installed and available."""
        return

    def _run(self):
        """Execute module."""
        # contact_executable = Path(FCC_path, self.params['executable'])

        contact_distance_cutoff = self.params["contact_distance_cutoff"]
        fraction_cutoff = self.params["fraction_cutoff"]
        threshold = self.params["threshold"]  # minimum number of models per cluster
        # strictness = self.params["strictness"]

        # Get the models generated in previous step
        models_to_cluster = self.previous_io.retrieve_models(individualize=True)

        log.info("Calculating contacts")
        contact_files = calculate_contacts(
            pdb_list=models_to_cluster, cutoff=contact_distance_cutoff, path=self.path  # type: ignore
        )

        log.info("Calculating the FCC matrix")
        matrix_outfile = Path("fcc.matrix")
        cluster_matrix = calculate_fcc_matrix(
            contact_list=contact_files, matrix_outfile=matrix_outfile
        )

        log.info("Clustering")
        cluster_check = False
        clusters = {}
        while not cluster_check:
            for threshold in range(self.params["threshold"], 0, -1):
                log.info(f"Clustering with threshold={threshold}")
                clusters = cluster_fcc(
                    matrix_fname=cluster_matrix,
                    similarity=fraction_cutoff,
                    minsize=threshold,
                    # cluster_outfile=cluster_outfile,
                )
                if not clusters:
                    log.info("[WARNING] No cluster was found, decreasing threshold!")
                else:
                    cluster_check = True
                    # pass the actual threshold back to the param dict
                    #  because it will be use in the detailed output
                    self.params["threshold"] = threshold
                    break
            if not cluster_check:
                # No cluster was obtained in any threshold
                cluster_check = True

        log.info("Writing clusters")

        # Prepare output and read the elements
        clt_dic = {}
        if clusters:
            cluster_info = collections.defaultdict(list)
            for k, v in clusters.items():
                cluster_info[v].append(k)
            # write the classic output file for compatibility reasons
            # log.info("Saving output to cluster.out")
            # cluster_out = Path("cluster.out")
            # with open(cluster_out, "w") as fh:
            #     cluster_fcc.output_clusters(fh, clusters)
            # fh.close()

            # clt_centers = {}
            for key in sorted(
                cluster_info, key=lambda k: len(cluster_info.get(k)), reverse=True
            ):
                cluster_id = key
                # cluster_center_id = clt.center.name - 1
                # cluster_center_pdb = models_to_cluster[cluster_center_id]

                clt_dic[cluster_id] = []
                # clt_centers[cluster_id] = cluster_center_pdb
                # clt_dic[cluster_id].append(cluster_center_pdb)

                for nb in cluster_info[key]:
                    model_id = nb + 1
                    model_pdb = models_to_cluster[model_id - 1]
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
            self.output_models = []
            for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
                cluster_id, _ = _e
                # sort the models by score
                clt_dic[cluster_id].sort()
                # rank the models
                for model_ranking, pdb in enumerate(clt_dic[cluster_id], start=1):
                    pdb.clt_id = cluster_id
                    pdb.clt_rank = cluster_rank
                    pdb.clt_model_rank = model_ranking
                    self.output_models.append(pdb)

            # Write unclustered structures
            write_structure_list(
                models_to_cluster, self.output_models, out_fname="clustfcc.tsv"
            )

            # Prepare clustfcc.txt
            output_fname = Path("clustfcc.txt")
            output_str = f"### clustfcc output ###{os.linesep}"
            output_str += os.linesep
            output_str += f"Clustering parameters {os.linesep}"
            output_str += (
                "> contact_distance_cutoff="
                f"{self.params['contact_distance_cutoff']}A"
                f"{os.linesep}"
            )
            output_str += (
                f"> fraction_cutoff={self.params['fraction_cutoff']}" f"{os.linesep}"
            )
            output_str += f"> threshold={self.params['threshold']}{os.linesep}"
            output_str += f"> strictness={self.params['strictness']}{os.linesep}"
            output_str += os.linesep
            output_str += (
                "Note: Models marked with * represent the center of the cluster"
                f"{os.linesep}"
            )
            output_str += f"-----------------------------------------------{os.linesep}"
            output_str += os.linesep
            output_str += f"Total # of clusters: {len(clusters)}{os.linesep}"

            for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
                cluster_id, _ = _e
                # center_pdb = clt_centers[cluster_id]
                model_score_l = [(e.score, e) for e in clt_dic[cluster_id]]
                model_score_l.sort()
                subset_score_l = [e[0] for e in model_score_l][:threshold]
                top_mean_score = np.mean(subset_score_l)
                top_std = np.std(subset_score_l)
                output_str += (
                    f"{os.linesep}"
                    "-----------------------------------------------"
                    f"{os.linesep}"
                    f"Cluster {cluster_rank} (#{cluster_id}, "
                    f"n={len(model_score_l)}, "
                    f"top{threshold}_avg_score = {top_mean_score:.2f} "
                    f"+-{top_std:.2f})"
                    f"{os.linesep}"
                )
                output_str += os.linesep
                output_str += f"clt_rank\tmodel_name\tscore{os.linesep}"
                for model_ranking, element in enumerate(model_score_l, start=1):
                    score, pdb = element
                    # if pdb.file_name == center_pdb.file_name:
                    #     output_str += (
                    #         f"{model_ranking}\t{pdb.file_name}\t{score:.2f}\t*"
                    #         f"{os.linesep}"
                    #     )
                    # else:
                    output_str += (
                        f"{model_ranking}\t{pdb.file_name}\t{score:.2f}" f"{os.linesep}"
                    )
            output_str += (
                "-----------------------------------------------" f"{os.linesep}"
            )

            log.info("Saving detailed output to clustfcc.txt")
            with open(output_fname, "w") as out_fh:
                out_fh.write(output_str)
        else:
            log.warning("No clusters were found")
            self.output_models = models_to_cluster

        self.output_models = models_to_cluster
        self.export_output_models()
