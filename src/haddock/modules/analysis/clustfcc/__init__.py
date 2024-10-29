"""Cluster modules with FCC.

The module takes the models generated in the previous step and calculates the
contacts between them. Then, the module calculates the FCC matrix and clusters
the models based on the calculated contacts.

For more details please check *Rodrigues, J. P. et al. Proteins: Struct. Funct. Bioinform. 80, 1810–1817 (2012)*
"""  # noqa: E501

import importlib.resources
import os
from pathlib import Path

from haddock import FCC_path, log
from haddock.core.defaults import CONTACT_FCC_EXEC, MODULE_DEFAULT_YAML
from haddock.core.typing import Union
from haddock.fcc import calc_fcc_matrix, cluster_fcc
from haddock.libs.libclust import (
    add_cluster_info,
    get_cluster_matrix_plot_clt_dt,
    plot_cluster_matrix,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libfcc import (
    calculate_pairwise_matrix,
    parse_contact_file,
    read_matrix,
    )
from haddock.libs.libsubprocess import JobInputFirst
from haddock.modules import BaseHaddockModule, get_engine, read_from_yaml_config
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.clustfcc.clustfcc import (
    get_cluster_centers,
    iterate_clustering,
    write_clusters,
    write_clustfcc_file,
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for clustering with FCC."""

    name = RECIPE_PATH.name

    def __init__(
        self,
        order: int,
        path: Path,
        initial_params: Union[Path, str] = DEFAULT_CONFIG,
    ) -> None:
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if FCC is installed and available."""
        # The FCC binary can be either in the default binary path or in the

        dcfg = read_from_yaml_config(DEFAULT_CONFIG)
        dcfg["executable"] = CONTACT_FCC_EXEC

    def _run(self) -> None:
        """Execute module."""
        contact_executable = Path(FCC_path, self.params["executable"])

        # Get the models generated in previous step
        models_to_clust = self.previous_io.retrieve_models(individualize=True)

        # Calculate the contacts for each model
        log.info("Calculating contacts")
        contact_jobs: list[JobInputFirst] = []
        for model in models_to_clust:
            pdb_f = Path(model.rel_path)  # type: ignore
            contact_f = Path(model.file_name.replace(".pdb", ".con"))  # type: ignore  # noqa : E501
            job = JobInputFirst(
                pdb_f,
                contact_f,
                CONTACT_FCC_EXEC,
                self.params["contact_distance_cutoff"],
            )
            contact_jobs.append(job)

        exec_mode = get_analysis_exec_mode(self.params["mode"])

        Engine = get_engine(exec_mode, self.params)
        engine = Engine(contact_jobs)
        engine.run()

        contact_file_l: list[str] = []
        not_found: list[str] = []
        for job in contact_jobs:
            if not job.output.exists():
                # NOTE: If there is no output, most likely the models are not in
                # contact there is no way of knowing how many models are not in
                # contact, it can be only one, or could be all of them.
                not_found.append(job.input.name)
                log.warning(f"Contact was not calculated for {job.input.name}")
            else:
                contact_file_l.append(str(job.output))

        if not_found:
            # No contacts were calculated, we cannot cluster
            self.finish_with_error("Several files were not generated:" f" {not_found}")

        log.info("Calculating the FCC matrix")
        parsed_contacts = parse_contact_file(
            contact_file_l,
            False,
        )

        # Imporant: matrix is a generator object, be careful with it
        matrix = calculate_pairwise_matrix(
            parsed_contacts,
            False,
        )

        # write the matrix to a file, so we can read it afterwards and don't
        #  need to reinvent the wheel handling this
        fcc_matrix_f = Path("fcc.matrix")
        with open(fcc_matrix_f, "w") as fh:
            for data in list(matrix):
                data_str = f"{data[0]} {data[1]} {data[2]:.2f} {data[3]:.3f}"
                data_str += os.linesep
                fh.write(data_str)

        # Cluster
        log.info("Clustering...")
        pool = read_matrix(
            fcc_matrix_f,
            self.params["clust_cutoff"],
            self.params["strictness"],
        )

        # iterate clustering until at least one cluster is found
        clusters, min_population = iterate_clustering(
            pool,
            self.params["min_population"],
        )
        self.params["min_population"] = min_population

        # Prepare output and read the elements
        if clusters:
            # Write the clusters
            write_clusters(clusters)

            # Get the cluster centers
            clt_dic, clt_centers = get_cluster_centers(
                clusters,
                models_to_clust,
            )

            # ranking clusters
            _scores, sorted_score_dic = rank_clusters(clt_dic, min_population)

            # Add this info to the models
            self.output_models = add_cluster_info(sorted_score_dic, clt_dic)

            # Write unclustered structures
            write_structure_list(
                models_to_clust,
                self.output_models,
                out_fname="clustfcc.tsv",
            )

            write_clustfcc_file(
                clusters, clt_centers, clt_dic, self.params, sorted_score_dic
            )
        else:
            log.warning("No clusters were found")
            self.output_models = models_to_clust  # type: ignore

        # Draw the matrix
        if self.params["plot_matrix"]:
            # Obtain final models indices
            final_order_idx, labels, cluster_ids = [], [], []
            for pdb in self.output_models:
                final_order_idx.append(models_to_clust.index(pdb))
                labels.append(pdb.file_name.replace(".pdb", ""))
                cluster_ids.append(pdb.clt_id)
            # Get custom cluster data
            matrix_cluster_dt, cluster_limits = get_cluster_matrix_plot_clt_dt(
                cluster_ids
            )

            # Define output filename
            html_matrix_basepath = "fcc_matrix"
            # Plot matrix
            html_matrixpath = plot_cluster_matrix(
                fcc_matrix_f,
                final_order_idx,
                labels,
                dttype="FCC",
                diag_fill=1,
                output_fname=html_matrix_basepath,
                matrix_cluster_dt=matrix_cluster_dt,
                cluster_limits=cluster_limits,
            )
            if html_matrixpath:
                log.info(f"Plotting matrix in {html_matrixpath}")
            else:
                log.warning("Cluster matrix was not generated")

        # Export models for next module
        self.export_io_models()
