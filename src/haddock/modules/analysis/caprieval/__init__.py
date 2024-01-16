"""Calculate CAPRI metrics."""
import json
from pathlib import Path

from haddock.gear.config import load as read_config
from haddock.core.typing import Any, FilePath
from haddock.modules import BaseHaddockModule, get_module_steps_folders
from haddock.modules import get_engine
from haddock.modules.analysis import get_analysis_exec_mode

from haddock.modules.analysis.caprieval.capri import (
    CAPRI,
    capri_cluster_analysis,
    merge_data,
    rearrange_ss_capri_output,
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")

CNS_MODULES = ["rigidbody",
               "flexref",
               "emscoring",
               "mdscoring",
               "mdref",
               "emref"]
WEIGHTS = ["w_elec", "w_vdw", "w_desolv", "w_bsa", "w_air"]


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

        # get the previous CNS step
        mod = len(sel_steps) - 2
        while mod > -1:
            st_name = sel_steps[mod].split("_")[1]
            if st_name in CNS_MODULES:
                cns_step = sel_steps[mod]
                break
            mod -= 1
        
        cns_params = read_config(Path("..", cns_step, "params.cfg"))
        key = list(cns_params['final_cfg'].keys())[0]
        scoring_pars = {kv: cns_params['final_cfg'][key][kv] for kv in WEIGHTS}
        
        scoring_params_fname = Path("weights_params.json")
        # write json file
        with open(scoring_params_fname, 'w', encoding='utf-8') as jsonf:
            json.dump(
                scoring_pars,
                jsonf,
                indent=4,
                )
        self.log(f"scoring parameters written in: {scoring_params_fname}")

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
        #  we can handle scenarios in wich the models are hetergoneous
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
