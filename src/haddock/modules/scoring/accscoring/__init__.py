"""Accessibility scoring module.

This module performs a solvent accessibility analysis based on some 
user-defined residues that should be buried or accessible.

If a supposedly buried (resp. accessible) residue is accessible (resp. buried),
the score should increase by one. The lower the final score the more consistent
the model with the user data.
"""
from pathlib import Path

from haddock.core.typing import FilePath
from haddock.modules import get_engine
from haddock.modules import BaseHaddockModule
from haddock.libs.libutil import parse_ncores
from haddock.libs.libparallel import get_index_list
from haddock.modules.scoring.accscoring.accscoring import (
    AccScore,
    AccScoreJob,
    prettify_df,
    rearrange_output,
)
from haddock.modules.analysis import (
    get_analysis_exec_mode,
    )

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to perform accessibility scoring."""

    name = RECIPE_PATH.name

    def __init__(self,
                 order: int,
                 path: Path,
                 initial_params: FilePath = DEFAULT_CONFIG) -> None:
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        return

    def _run(self) -> None:
        """Execute module."""
        try:
            models_to_score = self.previous_io.retrieve_models(
                individualize=True
                )
        except Exception as e:
            self.finish_with_error(e)

        nmodels = len(models_to_score)
        # index_list for the jobs with linear scaling
        ncores = parse_ncores(n=self.params['ncores'], njobs=nmodels)
        index_list = get_index_list(nmodels, ncores)

        # loading buried and accessible residue dictionaries
        buried_resdic = {
            key[-1]: value for key, value
            in self.params.items()
            if key.startswith("resdic_buried")
            }
        acc_resdic = {
            key[-1]: value for key, value
            in self.params.items()
            if key.startswith("resdic_accessible")
            }
        # remove _
        buried_resdic.pop("_")
        acc_resdic.pop("_")
        # finding the chains
        buried_violations_chains = [f"acc_{ch}" for ch in acc_resdic.keys()]
        acc_violations_chains = [f"bur_{ch}" for ch in buried_resdic.keys()]
        violations_chains = buried_violations_chains + acc_violations_chains

        # initialize jobs
        accscoring_jobs: list[AccScoreJob] = []
        
        for core in range(ncores):
            output_name = Path("accscoring_" + str(core) + ".tsv")
            viol_output_name = Path("violations_" + str(core) + ".tsv")
            accscore_obj = AccScore(
                model_list=models_to_score[index_list[core]:index_list[core + 1]],
                output_name=output_name,
                core=core,
                path=Path("."),
                buried_resdic=buried_resdic,
                acc_resdic=acc_resdic,
                cutoff=self.params["cutoff"],
                viol_output_name=viol_output_name,
                )
            
            job = AccScoreJob(
                accscore_obj,
            )
            accscoring_jobs.append(job)
        
        # Run accscoring Jobs
        exec_mode = get_analysis_exec_mode(self.params["mode"])
        Engine = get_engine(exec_mode, self.params)
        engine = Engine(accscoring_jobs)
        engine.run()

        # rearrange output
        output_name = Path("accscoring.tsv")
        accscoring_df = rearrange_output(
            output_name,
            path=Path("."),
            ncores=ncores
            )
        prettify_df(output_name, columns=["structure", "original_name", "md5", "score"], sortby="score")
        violations_df = rearrange_output(
            "violations.tsv",
            path=Path("."),
            ncores=ncores
            )
        prettify_df("violations.tsv", columns=["structure"] + violations_chains)

        self.output_models = models_to_score
        self.export_io_models()
