"""Compute contacts between chains in complexes.

The ``[contactmap]`` module aims at generating heatmaps and chordcharts of
the contacts observed in the input complexes.

If complexes are clustered, the analysis of contacts will be performed
based on all structures from each cluster.

**Heatmaps** are describing distances between all pairs of residues (both
intramolecular and intermolecular) for single models, and, for clusters,
the data selected with ``cluster_heatmap_datatype`` (by default the
probability of the shortest distance being under ``shortest_dist_threshold``).

**Chordcharts** are describing only intermolecular contacts in circles,
connecting with *chords* the two residues that are contacting. Residues are
colored by class (apolar, polar, positive, negative, nucleotide,
carbohydrate). Ribbons are colored when both residues share a class, or for
negative-positive pairs, and grey otherwise. Residue classes are not detected
interaction types. For clusters, a contact is drawn when the cluster
*average* distance is under threshold.

The ``ca-ca-dist`` is computed between CA atoms for amino acids, C4' atoms
for nucleotides and C1 atoms for carbohydrates.

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/analysis.html#contactmap-module>`_
"""

from copy import deepcopy
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath, SupportsRunT
from haddock.modules import BaseHaddockModule, get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.contactmap.contmap import (
    ClusteredContactMap,
    ContactsMap,
    get_clusters_sets,
    make_contactmap_report,
    topX_models,
)

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to compute complexes contacts and generate heatmap."""

    name = RECIPE_PATH.name

    def __init__(
        self,
        order: int,
        path: Path,
        *ignore: Any,
        init_params: FilePath = DEFAULT_CONFIG,
        **everything: Any,
    ) -> None:
        """Initialize class."""
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if contact executable is compiled."""
        return

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        models = []
        try:
            models = self.previous_io.retrieve_models(individualize=True)
        except AttributeError as e:
            self.finish_with_error(e)

        # Obtain clusters
        clusters_sets = get_clusters_sets(models)

        # Initiate holder of all jobs to be run by the `Scheduler`
        contact_jobs: list[SupportsRunT] = []
        # Loop over clusters
        for (clust_id, clust_rank), clt_models in clusters_sets.items():
            # In case of unclustered models
            if clust_id is None:
                # Obtain subset of top models
                top_models = topX_models(clt_models, topX=self.params["topX"])

                # Create single model analysis params
                single_models_params = deepcopy(self.params)
                single_models_params["single_model_analysis"] = True

                # Loop over models to analyse
                for model in top_models:
                    # Set names
                    modelfname = Path(model.file_name).stem
                    jobname = f"Unclustered_{modelfname}"
                    # Create a contact map object
                    contmap_job = ContactsMap(
                        Path(model.rel_path),
                        Path(jobname),
                        single_models_params,
                    )
                    contact_jobs.append(contmap_job)

            # For clustered models
            else:
                # Handles case where clustered models were not scored before
                rank_id = "Unranked" if clust_rank is None else clust_rank
                # Building basename for the job
                name = f"cluster{rank_id}"
                # Create a contact map object
                contmap_job = ClusteredContactMap(
                    [Path(model.rel_path) for model in clt_models],
                    Path(name),
                    self.params,
                )
                contact_jobs.append(contmap_job)

        # Find execution engine
        exec_mode = get_analysis_exec_mode(self.params["mode"])
        Engine = get_engine(exec_mode, self.params)
        engine = Engine(contact_jobs)
        engine.run()

        # Generate report
        make_contactmap_report(contact_jobs, "ContactMapReport.html")

        # Send models to the next step, no operation is done on them
        self.output_models = models
        self.export_io_models()
