"""Compute contacts between chains in complexes.

The ``[contactmap]`` module aims at generating heatmaps and chordcharts of
the contacts observed in the input complexes.

If complexes are clustered, the analysis of contacts will be performed
based on all structures from each cluster.

**Heatmaps** are describing the probability of contacts (<5A) between two
residues (both intramolecular and intermolecular).

**Chordcharts** are describing only intermolecular contacts in circles,
connecting with *chords* the two residues that are contacting.

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/analysis.html#contactmap-module>`_
"""

from copy import deepcopy
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, Iterator, FilePath, SupportsRunT
from haddock.modules import BaseHaddockModule
from haddock.modules import get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.contactmap.contmap import (
    ContactsMap,
    ClusteredContactMap,
    get_clusters_sets,
    make_contactmap_report,
    topX_models,
)
from haddock.libs.libutil import (
    get_available_memory,
    get_necessary_memory,
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
        if isinstance(self.previous_io, Iterator):
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)
        models = []
        try:
            models = self.previous_io.retrieve_models(individualize=True)
        except AttributeError as e:
            self.finish_with_error(e)

        # === IMPORTANT ================================================================
        # This modules uses a NxN distance matrix, this means that the memory
        # requirement will increase quadratically and can fail with an out-of-memory
        # error. Changing this behaviour would require a total re-write of the module
        # as of 04-2026 so instead we apply the following workaround:
        #  - Check what is the total size of the models (size is faster than reading)
        #  - Guesstimate how many atoms in total it would have based on the size
        #  - Calculate the expected matrix size and its memory requirements
        #  - Get how much memory the current host system has
        #  - If the system has less memory than needed, fail graciously
        current_memory = get_available_memory()
        needed_memory = get_necessary_memory(models) * self.params["ncores"]
        if current_memory < needed_memory:
            self.log(
                msg=(
                    f"Not enough memory to execute `contactmap` "
                    f"(needs {needed_memory:.2f}Gb has {current_memory:.2f}Gb). "
                    "! Skipping this module !"
                ),
                level="warning",
            )
            self.output_models = models
            self.export_io_models()
            return

        # ==============================================================================

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
