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


from os import linesep
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath, Optional, Union
from haddock.libs.libontology import PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.caprieval.capri import (
    CAPRI,
    capri_cluster_analysis,
    dump_weights,
    extract_data_from_capri_class,
    )
from pdbtools.pdb_wc import run as pdb_wc
from pdbtools.pdb_selmodel import run as pdb_selmodel


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

    def validate_reference(self, reference: Path) -> Path:
        """Validate the reference file by returning only one model.

        Parameters
        ----------
        reference : Path
            Path to the input reference structure, possibly containing
            an ensemble.

        Returns
        -------
        reference or first_model_path : Path
            Path to the reference structure to be used downstream.
        """
        # Extremly complicated stuff to manage the gathering of the sys.stdout,
        # as the pdb_tools.pdb_wc is basically writing on it.
        import sys
        from io import TextIOWrapper, BytesIO
        # Memorize previous sys.stdout
        original_stdout = sys.stdout
        # setup the new stdout environment
        sys.stdout = TextIOWrapper(BytesIO(), sys.stdout.encoding)

        # Count number of models
        with open(reference, "r") as fh:
            pdb_wc(fh, "m")
        # Get output
        sys.stdout.seek(0)  # Jump to the start
        wc_return = sys.stdout.read()  # Read output
        # Restore original stdout
        sys.stdout.close()
        sys.stdout = original_stdout
        # Parse output
        for line in wc_return.split("\n"):
            if "No. models" in line:
                sline = line.strip().split()
                nb_models = int(sline[-1])
                break
        # Return reference as only one structure present
        if nb_models == 1:
            return reference

        self.log(
            f"Multiple structures ({nb_models}) found in reference file. "
            "Using the first structure as reference. "
            "Consider providing different ensemble elements as "
            "reference_fname in sequential [caprieval] modules."
            )
        # Write new reference
        first_model_path = Path("first_model_reference.pdb")
        with open(reference, "r") as rin, open(first_model_path, "w") as rout:
            # Gather only first model (as string)
            for line in pdb_selmodel(rin, [1]):
                rout.write(line)
        return first_model_path

    def get_reference(self, models: list[PDBFile]) -> Path:
        """Manage to obtain the reference structure to be used downstream.

        Parameters
        ----------
        models : list[PDBFile]
            List of input model to be evaluated, among which the best
            can serve as reference structure if none provided.

        Returns
        -------
        reference : Path
            Path to the reference structure to be used downstream.
        """
        if self.params["reference_fname"] is not None:
            _reference = Path(self.params["reference_fname"])
            reference = self.validate_reference(_reference)
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
            reference = best_model_fname
        return reference

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
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

        # Get reference file
        reference = self.get_reference(models)

        # Each model is a job; this is not the most efficient way
        #  but by assigning each model to an individual job
        #  we can handle scenarios in which the models are hetergoneous
        #  for example during CAPRI scoring
        jobs: list[CAPRI] = []
        for i, model_to_be_evaluated in enumerate(models, start=1):
            # `models_to_be_evaluated` cannot be a list,
            # `CAPRI` class is expecting a single model
            if isinstance(model_to_be_evaluated, list):
                raise ValueError(
                    "CAPRI module cannot handle a list "
                    "of `model_to_be_evaluated`"
                )
            jobs.append(
                CAPRI(
                    identificator=i,
                    model=model_to_be_evaluated,
                    path=Path("."),
                    reference=reference,
                    params=self.params,
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

        extract_data_from_capri_class(
            capri_objects=jobs,
            output_fname=Path(".", "capri_ss.tsv"),
            sort_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
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
