"""Compute contacts between chains in complexes.

The ``[contactmap]`` module aims at generating heatmaps and chordcharts of
the contacts observed in the input complexes.

If complexes are clustered, the analysis of contacts will be performed
based on all structures from each cluster.

**Heatmaps** are describing the probability of contacts (<5A) between two
residues (both intramolecular and intermolecular).

**Chordcharts** are describing only intermolecular contacts in circles,
connecting with *chords* the two residues that are contacting.
"""

from copy import deepcopy
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath, SupportsRunT
from haddock.modules import BaseHaddockModule
from haddock.modules import get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.contactmap.contmap import (
    ContactsMap,
    ContactsMapJob,
    ClusteredContactMap,
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
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)
        try:
            models = self.previous_io.retrieve_models(individualize=True)
        except AttributeError as e:
            self.finish_with_error(e)

        # Obtain clusters
        clusters_sets = get_clusters_sets(models)

        # Initiate holder of all jobs to be run by the `Scheduler`
        contact_jobs: list[SupportsRunT] = []
        # Loop over clusters
        for clustid, clt_models in clusters_sets.items():
            # In case of unclustered models
            if clustid is None:
                # Obtain subset of top models
                top_models = topX_models(clt_models, topX=self.params["topX"])

                # Create single model analysis params
                single_models_params = deepcopy(self.params)
                single_models_params["single_model_analysis"] = True

                # Loop over models to analyse
                for model in top_models:
                    modelfname = Path(model.file_name).stem
                    # Create a job object
                    contmap_job = ContactsMapJob(
                        Path(f"Unclustered_contmap_{modelfname}"),
                        single_models_params,
                        modelfname,
                        # Create a contact map object
                        ContactsMap(
                            Path(model.rel_path),
                            Path(f"Unclustered_contmap_{modelfname}"),
                            single_models_params,
                            ),
                        )
                    contact_jobs.append(contmap_job)

            # For clustered models
            else:
                # Create a job object
                contmap_job = ContactsMapJob(
                    Path(f"cluster{clustid}_contmap"),
                    self.params,
                    f"Cluster_{clustid}",
                    # Create a contact map object
                    ClusteredContactMap(
                        [Path(model.rel_path) for model in clt_models],
                        Path(f"cluster{clustid}_contmap"),
                        self.params,
                        ),
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
