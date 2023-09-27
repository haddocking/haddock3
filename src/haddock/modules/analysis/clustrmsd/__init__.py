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
from pathlib import Path

import numpy as np

from haddock import log
from haddock.libs.libclust import (
    add_cluster_info,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libontology import ModuleIO
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_clusters,
    get_dendrogram,
    iterate_threshold,
    read_matrix,
    write_clusters,
    write_clustrmsd_file,
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with RMSD."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

        self.matrix_json = self._load_previous_io("rmsd_matrix.json")

    @classmethod
    def confirm_installation(cls):
        """Confirm if contact executable is compiled."""
        return

    def _run(self):
        """Execute module."""
        # Get the models generated in previous step
        models = self.previous_io.retrieve_models()

        # Cluster
        rmsd_matrix = read_matrix(
            self.matrix_json.input[0]
            )
                
        # getting clusters_list
        dendrogram = get_dendrogram(rmsd_matrix, self.params["linkage"])
        # setting tolerance
        if np.isnan(self.params["tolerance"]):
            self.log("tolerance is not defined")
            if self.params['criterion'] == "maxclust":
                tol = max(len(models) // 4 + 1, 2)
            else:
                tol = np.mean(dendrogram[:, 2])
            self.log(f"Setting tolerance to {tol:.2f} for criterion {self.params['criterion']}")  # noqa: E501
        else:
            if self.params['criterion'] == "maxclust":
                tol = int(self.params["tolerance"])
            elif self.params['criterion'] == "distance":
                tol = float(self.params["tolerance"])
            else:
                raise Exception(
                    f"unknown criterion {self.params['criterion']}"
                    )
        log.info(f"tolerance {tol}")
        self.params["tolerance"] = tol
        cluster_arr = get_clusters(
            dendrogram,
            self.params["tolerance"],
            self.params["criterion"]
            )

        # when crit == distance, apply clustering threshold
        if self.params['criterion'] == "distance":
            cluster_arr, threshold = iterate_threshold(
                cluster_arr,
                self.params['threshold']
                )
            self.params['threshold'] = threshold

        # print clusters
        unq_clusters = np.unique(cluster_arr)  # contains -1 (unclustered)
        clusters = [c for c in unq_clusters if c != -1]
        log.info(f"clusters = {clusters}")
        
        out_filename = Path('cluster.out')
        clt_dic, cluster_centers = write_clusters(
            clusters,
            cluster_arr,
            models,
            rmsd_matrix,
            out_filename,
            centers=True
            )

        # ranking clusters
        score_dic, sorted_score_dic = rank_clusters(
            clt_dic,
            self.params['threshold']
            )

        self.output_models = add_cluster_info(sorted_score_dic, clt_dic)
        
        # Write unclustered structures
        write_structure_list(models,
                             self.output_models,
                             out_fname="clustrmsd.tsv")

        write_clustrmsd_file(
            clusters,
            clt_dic,
            cluster_centers,
            score_dic,
            sorted_score_dic,
            self.params
            )

        self.export_io_models()
        # sending matrix to next step of the workflow
        matrix_io = ModuleIO()
        matrix_io.add(self.matrix_json.input[0])
        matrix_io.save(filename="rmsd_matrix.json")
