"""Cluster modules with FCC."""
import os
from pathlib import Path

from fcc.scripts import calc_fcc_matrix, cluster_fcc
from haddock import FCC_path, log
from haddock.core.typing import Union
from haddock.libs.libclust import (
    add_cluster_info,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libsubprocess import JobInputFirst
from haddock.modules import get_engine
from haddock.modules import BaseHaddockModule, read_from_yaml_config
from haddock.modules.analysis.clustfcc.clustfcc import (
    get_cluster_centers,
    iterate_clustering,
    write_clusters,
    write_clustfcc_file,
    )
from haddock.modules.analysis import get_analysis_exec_mode

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


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
        dcfg = read_from_yaml_config(DEFAULT_CONFIG)
        exec_path = Path(FCC_path, dcfg["executable"])

        if not os.access(exec_path, mode=os.F_OK):
            raise Exception(f"Required {str(exec_path)} file does not exist.")

        if not os.access(exec_path, mode=os.X_OK):
            raise Exception(f"Required {str(exec_path)} file is not executable")

        return

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
            contact_f = Path(model.file_name.replace(".pdb", ".con"))
            job = JobInputFirst(
                pdb_f,
                contact_f,
                contact_executable,
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
            self.finish_with_error(
                "Several files were not generated:"
                f" {not_found}"
                )

        log.info("Calculating the FCC matrix")
        parsed_contacts = calc_fcc_matrix.parse_contact_file(
            contact_file_l, False
            )

        # Imporant: matrix is a generator object, be careful with it
        matrix = calc_fcc_matrix.calculate_pairwise_matrix(
            parsed_contacts, False
            )

        # write the matrix to a file, so we can read it afterwards and don't
        #  need to reinvent the wheel handling this
        fcc_matrix_f = Path("fcc.matrix")
        with open(fcc_matrix_f, "w") as fh:
            for data in list(matrix):
                data_str = f"{data[0]} {data[1]} {data[2]:.2f} {data[3]:.3f}"
                data_str += os.linesep
                fh.write(data_str)
        fh.close()

        # Cluster
        log.info("Clustering...")
        pool = cluster_fcc.read_matrix(
            fcc_matrix_f,
            self.params["clust_cutoff"],
            self.params["strictness"],
            )

        # iterate clustering until at least one cluster is found
        clusters, min_population = iterate_clustering(
            pool,
            self.params['min_population'],
            )
        self.params['min_population'] = min_population

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
                clusters,
                clt_centers,
                clt_dic,
                self.params,
                sorted_score_dic
                )
        else:
            log.warning("No clusters were found")
            self.output_models = models_to_clust  # type: ignore

        self.export_io_models()
