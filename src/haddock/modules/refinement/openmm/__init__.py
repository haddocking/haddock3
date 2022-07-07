"""OpenMM refinement module for HADDOCK3."""
import os
import subprocess
from pathlib import Path

from haddock.libs.libontology import PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.modules import BaseHaddockModule
from haddock.modules.refinement.openmm.openmm import OPENMM


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


def runSubprocess(command_to_run):
    """Run subprocess."""
    outputOfSubprocess = subprocess.run([command_to_run],
                                        shell=True,
                                        capture_output=True,
                                        encoding='utf-8')
    return outputOfSubprocess.stdout.strip()


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 OpenMM module."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)
    
    def create_directories(self):
        """Create the necessary directories and provides the paths."""
        cwd = os.getcwd()

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
            dir_path = os.path.join(cwd, dir)
            os.mkdir(dir_path)
            directory_dict[dir] = dir_path

        return directory_dict

    @classmethod
    def confirm_installation(cls):
        """Confirm installation of openmm and pdfixer."""
        checkOpenMM = runSubprocess("conda list openmm --json")
        checkPdbfixer = runSubprocess("conda list pdbfixer --json")

        if(checkOpenMM == '[]'):
            raise Exception('OpenMM is not installed in conda.')
        if(checkPdbfixer == '[]'):
            raise Exception('OpenMM pdbfixer is not installed in conda.')
        return

    def _run(self):
        """Execute module."""
        previous_models = self.previous_io.retrieve_models()[0]
        self.log(f"previous models = {previous_models}")
        models_to_export = []
        # create directories
        directory_dict = self.create_directories()
        
        openmm_jobs = []
        for i, model_to_be_simulated in enumerate(previous_models, start=1):
            self.log(f"pdb {model_to_be_simulated}")
            self.log(f"filename {model_to_be_simulated.file_name}")
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
        ncores = self.params['ncores']
        openmm_engine = Scheduler(openmm_jobs, ncores=ncores)
        openmm_engine.run()

        # export models
        for pdb in os.listdir(directory_dict["openmm_output"]):
            pdbPath = os.path.join(directory_dict["openmm_output"], pdb)
            pdbToExport = PDBFile(pdbPath)
            models_to_export.append(pdbToExport)

        self.output_models = models_to_export
        self.export_output_models()

        self.log('Completed OpenMM module run.')
        topoaa_reminder = ['If you want to continue the haddock3 workflow',
                           ' after the OpenMM module, the next module should',
                           ' be topoaa, to rebuild the molecular topologies.'
                           ]
        self.log("".join(topoaa_reminder))
