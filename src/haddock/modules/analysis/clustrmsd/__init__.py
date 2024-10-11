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
* `n_clusters`: number of desired clusters (if `criterion` is `maxclust`).
* `clust_cutoff`: value of distance that separates distinct clusters (if `criterion` is
  ``distance``)
* `min_population` : analogously to the `clustfcc` module, it is the minimum number
  of models that should be present in a cluster to consider it. If criterion is
  `maxclust`, the value is ignored.

This module passes the path to the RMSD matrix is to the next step of the
workflow through the `rmsd_matrix.json` file, thus allowing to execute several
`clustrmsd` modules (possibly with different parameters) on the same RMSD
matrix.

.. _scipy routines: https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
"""  # noqa: E501
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Union
from haddock.libs.libclust import (
    add_cluster_info,
    clustrmsd_tolerance_params,
    get_cluster_matrix_plot_clt_dt,
    plot_cluster_matrix,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libontology import ModuleIO
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_clusters,
    get_dendrogram,
    get_matrix_path,
    iterate_min_population,
    order_clusters,
    read_matrix,
    write_clusters,
    write_clustrmsd_file,
    )

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


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
        min_population = self.params["min_population"]

        # getting clusters_list
        dendrogram = get_dendrogram(rmsd_matrix, self.params["linkage"])

        # adjust the parameters
        tolerance_param_name, tolerance = clustrmsd_tolerance_params(
            self.params,
            )
        
        log.info(
            f"Clustering with {tolerance_param_name} = {tolerance}, "
            f"and criterion {self.params['criterion']}"
            )
        
        cluster_arr = get_clusters(
            dendrogram,
            tolerance,
            self.params["criterion"],
            )

        # when crit == distance, apply clustering min_population
        if self.params['criterion'] == "distance":
            cluster_arr, min_population = iterate_min_population(
                cluster_arr,
                self.params['min_population'],
                )
            self.params['min_population'] = min_population
        
        clusters, cluster_arr = order_clusters(cluster_arr)
        log.info(f"clusters = {clusters}")
        
        out_filename = Path('cluster.out')
        clt_dic, cluster_centers = write_clusters(
            clusters,
            cluster_arr,
            models,
            rmsd_matrix,
            out_filename,
            centers=True,
            )

        # ranking clusters
        score_dic, sorted_score_dic = rank_clusters(
            clt_dic,
            self.params['min_population'],
            )

        self.output_models = add_cluster_info(sorted_score_dic, clt_dic)
        
        # Write unclustered structures
        write_structure_list(
            models,
            self.output_models,
            out_fname="clustrmsd.tsv",
            )  # type: ignore
        
        write_clustrmsd_file(
            clusters,
            clt_dic,
            cluster_centers,
            score_dic,
            sorted_score_dic,
            self.params,
            )

        # Draw the matrix
        if self.params['plot_matrix']:
            # Obtain final models indices
            final_order_idx, labels, cluster_ids = [], [], []
            for pdb in self.output_models:
                final_order_idx.append(models.index(pdb))
                labels.append(pdb.file_name.replace('.pdb', ''))
                cluster_ids.append(pdb.clt_id)
            # Get custom cluster data
            matrix_cluster_dt, cluster_limits = get_cluster_matrix_plot_clt_dt(
                cluster_ids
                )

            # Define output filename
            html_matrix_basepath = 'rmsd_matrix'
            # Plot matrix
            html_matrixpath = plot_cluster_matrix(
                get_matrix_path(self.matrix_json.input[0]),
                final_order_idx,
                labels,
                dttype='RMSD(Å)',
                reverse=True,
                diag_fill=0,
                output_fname=html_matrix_basepath,
                matrix_cluster_dt=matrix_cluster_dt,
                cluster_limits=cluster_limits,
                )
            if html_matrixpath:
                log.info(f"Plotting matrix in {html_matrixpath}")
            else:
                log.warning("Cluster matrix was not generated")

        self.export_io_models()
        # sending matrix to next step of the workflow
        matrix_io = ModuleIO()
        matrix_io.add(self.matrix_json.input[0])
        matrix_io.save(filename="rmsd_matrix.json")
