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
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath, Union
from haddock.libs.libontology import PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.modules import (
    BaseHaddockModule,
    get_engine,
    get_module_steps_folders,
    )
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.caprieval.capri import (
    CAPRI,
    capri_cluster_analysis,
    dump_weights,
    extract_data_from_capri_class,
    merge_data,
    rearrange_ss_capri_output,
)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to calculate the CAPRI metrics."""

    name = RECIPE_PATH.name

    def __init__(
        self,
        order: int,
        path: Path,
        *ignore: Any,
        init_params: FilePath = DEFAULT_CONFIG,
        **everything: Any,
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

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        models = self.previous_io.retrieve_models(individualize=True)
        if self.is_nested(models):
            raise ValueError(
                "CAPRI module cannot be executed after modules that produce a nested list of models"
            )

        # dump previously used weights
        dump_weights(self.order)

        # Sort by score to find the "best"
        models.sort()
        best_model = models[0]
        assert isinstance(best_model, PDBFile), "Best model is not a PDBFile"
        best_model_fname = best_model.rel_path

        if self.params["reference_fname"]:
            reference = Path(self.params["reference_fname"])
        else:
            self.log(
                "No reference was given. "
                "Using the structure with the lowest score from previous step"
            )
            reference = best_model_fname

        exec_mode = get_analysis_exec_mode(self.params["mode"])
        Engine = get_engine(exec_mode, self.params)

        _less_io = self.params["mode"] == "local" and not self.params["debug"]

        # Each model is a job; this is not the most efficient way
        #  but by assigning each model to an individual job
        #  we can handle scenarios in which the models are hetergoneous
        #  for example during CAPRI scoring
        jobs: list[CAPRI] = []
        for i, model_to_be_evaluated in enumerate(models, start=1):
            if isinstance(
                model_to_be_evaluated, list
            ):  # `models_to_be_evaluated` cannot be a list, `CAPRI` class is expecting a single model
                raise ValueError(
                    "CAPRI module cannot handle a list of `model_to_be_evaluated`"
                )
            jobs.append(
                CAPRI(
                    identificator=str(i),
                    model=model_to_be_evaluated,
                    path=Path("."),
                    reference=reference,
                    params=self.params,
                    debug=not _less_io,
                )
            )

        engine = Engine(jobs)
        engine.run()

        if _less_io and isinstance(engine, Scheduler):
            jobs = engine.results
            extract_data_from_capri_class(
                capri_objects=jobs,
                output_fname=Path(".", "capri_ss.tsv"),
                sort_key=self.params["sortby"],
                sort_ascending=self.params["sort_ascending"],
            )

        else:
            self.log(
                msg=(
                    "DEPRECATION NOTICE: This execution mode (debug=True) "
                    "will no longer be supported in the next version."
                    ),
                level="warning",
            )
            jobs = merge_data(jobs)

            # Each job created one .tsv, unify them:
            rearrange_ss_capri_output(
                output_name="capri_ss.tsv",
                output_count=len(jobs),
                sort_key=self.params["sortby"],
                sort_ascending=self.params["sort_ascending"],
                path=Path("."),
            )

        capri_cluster_analysis(
            capri_list=jobs,
            model_list=models,  # type: ignore # ignore this here only if we are checking the return type of `retrieve_models` is not nested!!
            output_fname="capri_clt.tsv",
            clt_threshold=self.params["clt_threshold"],
            # output_count=len(capri_jobs),
            sort_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            path=Path("."),
        )

        # Send models to the next step,
        #  no operation is done on them
        self.output_models = models  # type: ignore # ignore this here only if we are checking the return type of `retrieve_models` is not nested!!
        self.export_io_models()
