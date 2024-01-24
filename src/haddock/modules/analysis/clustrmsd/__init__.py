"""
RMSD clustering module.

This module takes in input the RMSD matrix calculated in the previous step and
performs a hierarchical clustering procedure on it, leveraging `scipy routines`_
for this purpose.

Essentially, the procedure amounts at lumping the input models in a
progressively coarser hierarchy of clusters, called the dendrogram.

Four parameters can be defined in this context:

* `linkage`: governs the way clusters are merged together in the creation of
  the dendrogram
* `criterion`: defines the prescription to cut the dendrogram and obtain the
  desired clusters
* `tolerance`:
    * if `criterion` is ``maxclust``, this is the number of desired
      clusters.
    * if `criterion` is ``distance``, it must be the value of
      distance that separates distinct clusters.

  If not specified, the default is calculated either as the total number of
  models divided by four (`maxclust`) or as the average of the dendrogram height
  (`distance`)
* `threshold` : analogously to the `clustfcc` module, it is the minimum number
  of models that should be present in a cluster to consider it

This module passes the path to the RMSD matrix is to the next step of the
workflow through the `rmsd_matrix.json` file, thus allowing to execute several
`clustrmsd` modules (possibly with different parameters) on the same RMSD
matrix.

.. _scipy routines: https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
"""  # noqa: E501
import os
from pathlib import Path

import numpy as np

from haddock import log

# from haddock.core.typing import FilePath
from haddock.libs.libclust import write_structure_list, plot_cluster_matrix
from haddock.libs.libontology import ModuleIO
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_cluster_center,
    get_clusters,
    get_dendrogram,
    get_matrix_path,
    iterate_threshold,
    read_matrix,
    )
from typing import Union


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with RMSD."""

    name = RECIPE_PATH.name

    def __init__(
            self,
            order: int,
            path: Path,
            initial_params: Union[Path, str] = DEFAULT_CONFIG,
            ) -> None:
        super().__init__(order, path, initial_params)

        self.matrix_json = self._load_previous_io("rmsd_matrix.json")

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if contact executable is compiled."""
        return

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        models = self.previous_io.retrieve_models()
        # Cluster
        rmsd_matrix = read_matrix(self.matrix_json.input[0])
        # loading parameters
        linkage_type = self.params["linkage"]
        crit = self.params["criterion"]
        threshold = self.params["threshold"]

        # getting clusters_list
        dendrogram = get_dendrogram(rmsd_matrix, linkage_type)
        # setting tolerance
        tol: Union[float, int]
        if np.isnan(self.params["tolerance"]):
            self.log("tolerance is not defined")
            if crit == "maxclust":
                tol = max(len(models) // 4 + 1, 2)
            else:
                tol = np.mean(dendrogram[:, 2])
            self.log(f"Setting tolerance to {tol:.2f} for criterion {crit}")
        else:
            if crit == "maxclust":
                tol = int(self.params["tolerance"])
            elif crit == "distance":
                tol = float(self.params["tolerance"])
            else:
                raise Exception(f"unknown criterion {crit}")
        log.info(f"tolerance {tol}")
        cluster_arr = get_clusters(dendrogram, tol, crit)

        # when crit == distance, apply clustering threshold
        if crit == "distance":
            cluster_arr = iterate_threshold(cluster_arr, threshold)

        # print clusters
        unq_clusters = np.unique(cluster_arr)  # contains -1 (unclustered)
        clusters = [c for c in unq_clusters if c != -1]
        log.info(f"clusters = {clusters}")

        # initialising cluster centers
        n_obs = len(cluster_arr)
        cluster_centers = {}
        # preparing output
        clt_dic = {}
        log.info("Saving output to cluster.out")
        cluster_out = Path("cluster.out")
        with open(cluster_out, "w") as fh:
            for cl_id in clusters:
                if cl_id != -1:
                    npw = np.where(cluster_arr == cl_id)[0]
                    clt_dic[cl_id] = [models[n] for n in npw]
                    fh.write(f"Cluster {cl_id} -> ")
                    clt_center = get_cluster_center(npw, n_obs, rmsd_matrix)
                    cluster_centers[cl_id] = models[clt_center].file_name  # type: ignore  # noqa : E501
                    for el in npw[:-1]:
                        fh.write(f"{el + 1} ")
                    fh.write(f"{npw[-1] + 1}")
                    fh.write(os.linesep)
        # rank the clusters
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
                pdb.clt_id = int(cluster_id)
                pdb.clt_rank = cluster_rank
                pdb.clt_model_rank = model_ranking
                self.output_models.append(pdb)

        # Write unclustered structures
        write_structure_list(models, self.output_models, out_fname="clustrmsd.tsv")  # type: ignore  # noqa : E501

        # Prepare clustrmsd.txt
        output_fname = Path("clustrmsd.txt")
        output_str = f"### clustrmsd output ###{os.linesep}"
        output_str += os.linesep
        output_str += f"Clustering parameters {os.linesep}"
        output_str += f"> linkage_type={linkage_type}{os.linesep}"
        output_str += f"> criterion={crit}{os.linesep}"
        output_str += f"> tolerance={tol:.2f}{os.linesep}"
        output_str += f"> threshold={threshold}{os.linesep}"
        output_str += os.linesep

        output_str += f"{'-' * 47}{os.linesep}"
        output_str += os.linesep
        output_str += f"Total # of clusters: {len(clusters)}{os.linesep}"
        for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
            cluster_id, _ = _e

            model_score_l = [(e.score, e) for e in clt_dic[cluster_id]]
            model_score_l.sort()
            top_score = score_dic[cluster_id]

            output_str += (
                f"{os.linesep}"
                f"{'-' * 47}"
                f"{os.linesep}"
                f"Cluster {cluster_rank} (#{cluster_id}, "
                f"n={len(model_score_l)}, "
                f"top{threshold}_avg_score = {top_score:.2f})"
                f"{os.linesep}"
                )
            output_str += os.linesep
            output_str += f"clt_rank\tmodel_name\tscore{os.linesep}"
            for model_ranking, element in enumerate(model_score_l, start=1):
                score, pdb = element
                # is the model the cluster center?
                if pdb.file_name == cluster_centers[cluster_id]:
                    output_str += (
                        f"{model_ranking}\t{pdb.file_name}\t{score:.2f}\t*"
                        f"{os.linesep}"
                        )
                else:
                    output_str += (
                        f"{model_ranking}\t{pdb.file_name}\t{score:.2f}"
                        f"{os.linesep}"
                        )
        output_str += f"{'-' * 47}{os.linesep}"
        log.info("Saving detailed output to clustrmsd.txt")
        with open(output_fname, "w") as out_fh:
            out_fh.write(output_str)

        # Draw the matrix
        if self.params['plot_matrix']:
            # Obtain final models indices
            final_order_idx, labels = [], []
            for pdb in self.output_models:
                final_order_idx.append(models.index(pdb))
                labels.append(pdb.file_name.replace('.pdb', ''))

            # Define output filename
            html_matrix_basepath = 'rmsd_matrix'
            # Plot matrix
            html_matrixpath = plot_cluster_matrix(
                get_matrix_path(self.matrix_json.input[0]),
                final_order_idx,
                labels,
                dttype='RMSD (Å)',
                reverse=True,
                diag_fill=0,
                output_fname=html_matrix_basepath,
                )
            log.info(f"Plotting matrix in {html_matrixpath}")

        self.export_io_models()
        # sending matrix to next step of the workflow
        matrix_io = ModuleIO()
        matrix_io.add(self.matrix_json.input[0])
        matrix_io.save(filename="rmsd_matrix.json")
