"""
Surface accessibility scoring module.

This module performs a solvent accessibility analysis based on 
some user-defined residues that should be buried or accessible.

If a supposedly buried (resp. accessible) residue is accessible (resp. buried),
the score should increase by one. The lower the final score the more consistent
the model with the user data.

To run this module, the user must provide a dictionary with the buried and/or accessible
residues. The keys of the dictionary should be the chain identifiers and the values
should be lists of residue numbers.

Example:

>>> resdic_buried_A: [1, 2, 3]
>>> resdic_accessible_B: [4, 5, 6]
"""
from pathlib import Path
import os

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath
from haddock.modules import get_engine
from haddock.modules.scoring import ScoringModule
from haddock.libs.libontology import PDBFile
from haddock.libs.libutil import parse_ncores
from haddock.libs.libparallel import get_index_list, Scheduler
from haddock.modules.scoring.sasascore.sasascore import (
    AccScore,
    extract_data_from_accscore_class
    )

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(ScoringModule):
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
        if self.params["mode"] == "mpi":
            input_ncores_slurm = int(os.getenv('SLURM_CPUS_ON_NODE', 1))
            input_ncores_pbs = int(os.getenv('PBS_CPUS_ON_NODE', 1))
            ncores = max(input_ncores_slurm, input_ncores_pbs)
        else:
            ncores = self.params["ncores"]
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
        
        self.output_models: list[PDBFile] = []
        # initialize jobs
        sasascore_jobs: list[AccScore] = []
        for i, model_to_be_evaluated in enumerate(models_to_score):
            accscore_obj = AccScore(
                model=model_to_be_evaluated,
                path=Path("."),
                buried_resdic=buried_resdic,
                acc_resdic=acc_resdic,
                cutoff=self.params["cutoff"],
                probe_radius=self.params["probe_radius"],
                identificator=i,
                )
            sasascore_jobs.append(accscore_obj)
            # append model to output models
            self.output_models.append(model_to_be_evaluated)
        
        # Run sasascore Jobs using Scheduler
        engine = Scheduler(
            ncores=ncores,
            tasks=sasascore_jobs,
        )
        engine.run()

        # extract results and overwrite scores
        sasascore_jobs = engine.results
        sasascore_jobs = sorted(sasascore_jobs,
                                key=lambda accscore: accscore.identificator)
        
        for i, job in enumerate(sasascore_jobs):
            pdb = self.output_models[i]
            pdb.score = job.data[3]
        output_name = Path("sasascore.tsv")
        self.output(output_name)

        # now violations
        viol_output_name = Path("violations.tsv")
        extract_data_from_accscore_class(
            sasascore_objects=sasascore_jobs,
            violations_output_fname=viol_output_name
        )

        self.export_io_models()
