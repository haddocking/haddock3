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
from haddock.libs.libontology import ModuleIO
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_clusters,
    get_dendrogram,
    read_matrix,
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
        linkage_type = self.params["linkage"]
        dendrogram = get_dendrogram(rmsd_matrix, linkage_type)
        crit = self.params["criterion"]
        if np.isnan(self.params["tolerance"]):
            self.log("tolerance is not defined")
            if crit == "maxclust":
                tol = len(models) // 4 + 1
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
        cluster_list = get_clusters(dendrogram, tol, crit)
        clusters = np.unique(cluster_list)
        log.info(f"clusters = {clusters}")
        # TODO: getting cluster centers. is this really necessary?

        # preparing output
        clt_dic = {}
        log.info('Saving output to cluster.out')
        cluster_out = Path('cluster.out')
        with open(cluster_out, 'w') as fh:
            for cl_id in clusters:
                fh.write(f"Cluster {cl_id} -> ")
                npw = np.where(cluster_list == cl_id)
                clt_dic[cl_id] = [models[n] for n in npw[0]]
                for el in npw[0][:-1]:
                    fh.write(str(el + 1) + " ")
                fh.write(str(npw[0][-1] + 1))
                fh.write(os.linesep)
        # rank the clusters
        threshold = self.params['threshold']
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
            for model_ranking, pdb in enumerate(clt_dic[cluster_id],
                                                start=1):
                pdb.clt_id = cluster_id
                pdb.clt_rank = cluster_rank
                pdb.clt_model_rank = model_ranking
                self.output_models.append(pdb)
        
        # Prepare clustrmsd.txt
        output_fname = Path('clustrmsd.txt')
        output_str = f'### clustrmsd output ###{os.linesep}'
        output_str += os.linesep
        output_str += f'Clustering parameters {os.linesep}'
        output_str += f"> linkage_type={linkage_type}{os.linesep}"
        output_str += f"> criterion={crit}{os.linesep}"
        output_str += f"> tolerance={tol:.2f}{os.linesep}"
        output_str += f"> threshold={threshold}{os.linesep}"
        output_str += os.linesep
        
        output_str += (
            f"-----------------------------------------------{os.linesep}")
        output_str += os.linesep
        output_str += f'Total # of clusters: {len(clusters)}{os.linesep}'
        for cluster_rank, _e in enumerate(sorted_score_dic, start=1):
            cluster_id, _ = _e
            
            model_score_l = [(e.score, e) for e in clt_dic[cluster_id]]
            model_score_l.sort()
            top_score = score_dic[cluster_id]

            output_str += (
                f"{os.linesep}"
                "-----------------------------------------------"
                f"{os.linesep}"
                f"Cluster {cluster_rank} (#{cluster_id}, "
                f"n={len(model_score_l)}, "
                f"top{threshold}_avg_score = {top_score:.2f})"
                f"{os.linesep}")
            output_str += os.linesep
            output_str += f'clt_rank\tmodel_name\tscore{os.linesep}'
            for model_ranking, element in enumerate(model_score_l, start=1):
                score, pdb = element
                
                output_str += (
                    f"{model_ranking}\t{pdb.file_name}\t{score:.2f}"
                    f"{os.linesep}")
        output_str += (
            "-----------------------------------------------"
            f"{os.linesep}")
        log.info('Saving detailed output to clustrmsd.txt')
        with open(output_fname, 'w') as out_fh:
            out_fh.write(output_str)

        self.export_output_models()
        # sending matrix to next step of the workflow
        matrix_io = ModuleIO()
        matrix_io.add(self.matrix_json.input[0])
        matrix_io.save(filename="rmsd_matrix.json")
