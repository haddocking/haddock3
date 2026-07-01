"""Calculate CAPRI metrics for the input models.

By default the following metrics are calculated:

- FNAT (fraction of native contacts), namely the fraction of
    intermolecular contacts in the docked complex that are also
    present in the reference complex.
- IRMSD (interface root mean square deviation), namely the RMSD
    of the interface of the docked complex with respect
    to the reference complex.
- LRMSD (ligand root mean square deviation), namely the RMSD of the
    ligand of the docked complex with respect to the
    reference complex upon superposition of the receptor.
- DOCKQ, a measure of the quality of the docked model obtained
    by combining FNAT, IRMSD and LRMSD (see
    Basu and Wallner 2016,  11 (8), e0161879).
- ILRMSD (interface ligand root mean square deviation), the RMSD of the
    ligand of the docked complex with respect to the reference complex
    upon superposition of the interface of the receptor.
- GLOBAL_RMSD, the full RMSD between the reference and the model.

The following files are generated:

- **capri_ss.tsv**: a table with the CAPRI metrics for each model.
- **capri_clt.tsv**: a table with the CAPRI metrics for each cluster of models (if clustering information is available).

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/analysis.html#caprieval-module>`_
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath, Iterator, Union
from haddock.libs.libontology import PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libaa2cg import martinize
from haddock.libs.libstructure import find_ff
from haddock.libs.libpdb import handle_input_reference
from haddock.modules import BaseHaddockModule
from haddock.libs.libcapri import (
    CAPRI,
    capri_cluster_analysis,
    extract_data_from_capri_class,
    extract_models_best_references,
)
from haddock.modules.analysis.caprieval.capri import dump_weights


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to calculate the CAPRI metrics."""

    name = RECIPE_PATH.name

    def __init__(
        self,
        order: int,
        path: Path,
        init_params: FilePath = DEFAULT_CONFIG,
    ) -> None:
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if contact executable is compiled."""
        return

    @staticmethod
    def is_nested(models: list[Union[PDBFile, list[PDBFile]]]) -> bool:
        for model in models:
            if isinstance(model, list):
                return True
        return False

    def get_reference(self, models: list[PDBFile]) -> list[Path]:
        """Manage to obtain the reference structure to be used downstream.

        Parameters
        ----------
        models : list[PDBFile]
            List of input model to be evaluated, among which the best
            can serve as reference structure if none provided.

        Returns
        -------
        references : list[Path]
            List of paths to the reference(s) structure to be used downstream.
        """
        if self.params["reference_fname"]:
            _reference = Path(self.params["reference_fname"])
            references = handle_input_reference(_reference)
        else:
            self.log(
                "No reference structure provided. "
                "Using the structure with the lowest score from previous step"
            )
            # Sort by score to find the "best"
            models.sort()
            best_model = models[0]
            assert isinstance(best_model, PDBFile), "Best model is not a PDBFile"
            best_model_fname = best_model.rel_path
            references = [best_model_fname]
        return references

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        if isinstance(self.previous_io, Iterator):
            # if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        models = self.previous_io.retrieve_models(individualize=True)
        if self.is_nested(models):
            raise ValueError(
                "CAPRI module cannot be executed after "
                "modules that produce a nested list of models."
            )

        # dump previously used weights
        dump_weights(self.order)

        # Find force-field
        ff = find_ff(models)
        # Get reference file
        if ff == "martini2":
            references = [
                Path(martinize(ref_aa, self.path.resolve().parent, False))
                for ref_aa in self.get_reference(models)
            ]
        else:
            references = self.get_reference(models)

        # Each model is a job; this is not the most efficient way
        #  but by assigning each model to an individual job
        #  we can handle scenarios in which the models are hetergoneous
        #  for example during CAPRI scoring
        jobs: list[CAPRI] = []
        # Loop over models
        for i, model_to_be_evaluated in enumerate(models, start=1):
            # `models_to_be_evaluated` cannot be a list,
            # `CAPRI` class is expecting a single model
            if isinstance(model_to_be_evaluated, list):
                raise ValueError(
                    "CAPRI module cannot handle a list of `model_to_be_evaluated`"
                )
            # Loop over references
            for ref_id, reference in enumerate(references, start=1):
                jobs.append(
                    CAPRI(
                        identificator=i,
                        model=model_to_be_evaluated,
                        path=Path("."),
                        reference=reference,
                        params=self.params,
                        ref_id=ref_id,
                        ff=ff,
                    )
                )

        engine = Scheduler(
            tasks=jobs,
            ncores=self.params["ncores"],
            max_cpus=self.params["max_cpus"],
        )
        engine.run()

        jobs = engine.results
        jobs = sorted(jobs, key=lambda capri: capri.identificator)

        # Extract best references per input model
        best_ref_jobs = extract_models_best_references(jobs)

        # Write standard capri_ss file
        extract_data_from_capri_class(
            capri_objects=best_ref_jobs,
            output_fname=Path(".", "capri_ss.tsv"),
            sort_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            add_reference_id=len(references) > 1,
        )

        # Perform cluster analysis
        capri_cluster_analysis(
            capri_list=best_ref_jobs,
            model_list=models,  # type: ignore # ignore this here only if we are checking the return type of `retrieve_models` is not nested!!
            output_fname="capri_clt.tsv",
            clt_threshold=self.params["clt_threshold"],
            # output_count=len(capri_jobs),
            sort_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            path=Path("."),
        )

        # In case multiple references are provided, generate an additional file
        # containing the information to traceback metrics related to each ref
        if len(references) > 1:
            extract_data_from_capri_class(
                capri_objects=jobs,
                output_fname=Path(".", "capri_ss_multiref.tsv"),
                sort_key=self.params["sortby"],
                sort_ascending=self.params["sort_ascending"],
                add_reference_id=True,
            )

        # Send models to the next step, no operation is done on them
        self.output_models = models  # type: ignore # ignore this here only if we are checking the return type of `retrieve_models` is not nested!!
        self.export_io_models()
