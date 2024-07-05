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

The following files are generated:

- **capri_ss.tsv**: a table with the CAPRI metrics for each model.
- **capri_clt.tsv**: a table with the CAPRI metrics for each cluster of models (if clustering information is available).
"""
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.modules import BaseHaddockModule, get_module_steps_folders
from haddock.modules import get_engine
from haddock.modules.analysis import get_analysis_exec_mode

from haddock.modules.analysis.caprieval.capri import (
    CAPRI,
    capri_cluster_analysis,
    get_previous_cns_step,
    merge_data,
    rearrange_ss_capri_output,
    save_scoring_weights,
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to calculate the CAPRI metrics."""

    name = RECIPE_PATH.name

    def __init__(self,
                 order: int,
                 path: Path,
                 *ignore: Any,
                 init_params: FilePath = DEFAULT_CONFIG,
                 **everything: Any) -> None:
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

        models = self.previous_io.retrieve_models(
            individualize=True
            )
        
        # dump previously used weights
        sel_steps = get_module_steps_folders(Path(".."))
        cns_step = get_previous_cns_step(sel_steps, self.order)
        
        if cns_step:
            self.log(f"Found previous CNS step: {cns_step}")
            scoring_params_fname = save_scoring_weights(cns_step)
            self.log(f"Saved scoring weights to: {scoring_params_fname}")
        else:
            self.log("No previous CNS step found. Cannot save scoring weights.")

        # Sort by score to find the "best"
        models.sort()
        best_model_fname = Path(models[0].rel_path)

        if self.params["reference_fname"]:
            reference = Path(self.params["reference_fname"])
        else:
            self.log(
                "No reference was given. "
                "Using the structure with the lowest score from previous step")
            reference = best_model_fname

        # Each model is a job; this is not the most efficient way
        #  but by assigning each model to an individual job
        #  we can handle scenarios in which the models are hetergoneous
        #  for example during CAPRI scoring
        jobs: list[CAPRI] = []
        for i, model_to_be_evaluated in enumerate(models, start=1):
            jobs.append(
                CAPRI(
                    identificator=str(i),
                    model=model_to_be_evaluated,
                    path=Path("."),
                    reference=reference,
                    params=self.params
                    )
                )
        
        exec_mode = get_analysis_exec_mode(self.params["mode"])

        Engine = get_engine(exec_mode, self.params)
        engine = Engine(jobs)
        engine.run()

        # very ugly way of loading the capri metrics back into
        #  the CAPRI object, there's definitively a better way
        #  of doing this
        jobs = merge_data(jobs)

        # Each job created one .tsv, unify them:
        rearrange_ss_capri_output(
            output_name="capri_ss.tsv",
            output_count=len(jobs),
            sort_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            path=Path(".")
            )

        capri_cluster_analysis(
            capri_list=jobs,
            model_list=models,
            output_fname="capri_clt.tsv",
            clt_threshold=self.params["clt_threshold"],
            # output_count=len(capri_jobs),
            sort_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            path=Path(".")
            )

        # Send models to the next step,
        #  no operation is done on them
        self.output_models = models
        self.export_io_models()
