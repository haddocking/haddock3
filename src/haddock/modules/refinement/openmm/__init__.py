"""OpenMM refinement module for HADDOCK3."""
import os
from contextlib import suppress
from pathlib import Path

from haddock.libs.libontology import PDBFile
from haddock.libs.libsubprocess import run_subprocess
from haddock.modules import BaseHaddockModule, get_engine


# allow general testing when OpenMM is not installed
with suppress(ImportError):
    from haddock.modules.refinement.openmm.openmm import OPENMM


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 OpenMM module."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

    def create_directories(self):
        """Create the necessary directories and provides the paths."""
        directory_list = [
            "pdbfixer",
            "solvation_boxes",
            "intermediates",
            "md_raw_output",
            "openmm_output"
            ]

        directory_dict = {}
        for dir in directory_list:
            self.log(f"Creating directory {dir}")
            os.mkdir(dir)
            directory_dict[dir] = dir
        
        return directory_dict

    @classmethod
    def confirm_installation(cls):
        """Confirm installation of openmm and pdfixer."""
        checkOpenMM = run_subprocess("conda list openmm --json")
        checkPdbfixer = run_subprocess("conda list pdbfixer --json")

        if (checkOpenMM == '[]'):
            raise Exception('OpenMM is not installed in conda.')
        if (checkPdbfixer == '[]'):
            raise Exception('OpenMM pdbfixer is not installed in conda.')
        return

    def _run(self):
        """Execute module."""
        previous_models = self.previous_io.retrieve_models(
            individualize=True
            )

        previous_models.sort()
        models_to_export = []
        # create directories
        directory_dict = self.create_directories()

        openmm_jobs = []
        for i, model_to_be_simulated in enumerate(previous_models, start=1):
            openmm_jobs.append(
                OPENMM(
                    identificator=i,
                    model=model_to_be_simulated,
                    path=Path("."),
                    params=self.params,
                    directory_dict=directory_dict
                    )
                )

        # running jobs
        # currently it only accepts local and mpi.
        # And at the end of the execution mpi fails
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(openmm_jobs)
        engine.run()

        # export models
        output_pdbs = os.listdir(directory_dict["openmm_output"])

        for pdb in output_pdbs:
            pdbPath = os.path.join(directory_dict["openmm_output"], pdb)
            pdbToExport = PDBFile(pdbPath)
            models_to_export.append(pdbToExport)
        
        self.output_models = models_to_export
        self.export_output_models()

        self.log('Completed OpenMM module run.')
        self.log('If you want to continue the haddock3 workflow',
                 ' after the OpenMM module, the next module should',
                 ' be topoaa, to rebuild the molecular topologies.')
