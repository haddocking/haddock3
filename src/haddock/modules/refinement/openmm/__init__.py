"""OpenMM refinement module for HADDOCK3."""
from haddock.libs.libontology import Format, PDBFile
from haddock.modules import BaseHaddockModule
import os
import subprocess
from pathlib import Path
#import haddock.modules.refinement.openmm.openmmfunctions as openmmfunctions
# my additions
from haddock.modules.refinement.openmm.openmm import OPENMM
from haddock.libs.libparallel import Scheduler

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")

class HaddockModule(BaseHaddockModule):
    """HADDOCK3 OpenMM module."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)
    
    def create_directories(self):
        """Creates the necessary directories and provides the paths."""
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
        def runSubprocess(command_to_run):
            outputOfSubprocess = subprocess.run([command_to_run], shell=True, capture_output=True, encoding='utf-8')
            return outputOfSubprocess.stdout.strip()
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

        #openmmPdbfixer_output_Directory, modeller_solvationbox_pdbs_Directory, openmm_intermediate_structures, openmm_md_raw_output_Directory, openmm_output = self.createDirectories()
        directory_dict = self.create_directories()
        
        openmm_jobs = []
        for i, model_to_be_simulated in enumerate(previous_models, start=1):
            self.log(f"pdb {model_to_be_simulated}")
            #self.log(f"pdb {model_to_be_evaluated.keys()}")
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

        #for pdb in previous_models:
        #    if(pdb.file_type is Format.PDB):
        #        pdb_filePath = f'{Path(pdb.path) / pdb.file_name}'
        #        pdb_filepath_openmm_pdbfixer_directory = os.path.join(directory_dict["pdbfixer"], pdb.file_name)                
        #        self.log('Fixing pdb with openmm pdbfixer.')
        #        openmmfunctions.OpenmmPdbfixer(self, pdb_filePath, pdb.file_name, directory_dict["pdbfixer"])
        #        if(self.params['implicit_solvent'] == False):
        #            self.log(f'Building solvation box for file: {pdb.file_name}')
        #            contains_xray_crystallography_cell_data = openmmfunctions.Does_pdb_contain_xray_crystallography_cell_data(pdb_filePath)
        #            openmmfunctions.createSolvationBox(self, pdb_filepath_openmm_pdbfixer_directory, pdb.file_name, directory_dict["solvation_boxes"], self.params['forcefield'][0], self.params['explicit_solvent_model'][0], contains_xray_crystallography_cell_data)
#
        #self.log(f'Starting openMM simulations.')
        #if(self.params['implicit_solvent']):
        #    for pdb in os.listdir(directory_dict["pdbfixer"]):
        #        pdbPath = os.path.join(directory_dict["pdbfixer"], pdb)
        #        self.log(f'Starting openMM with file: {pdbPath}')
        #        openmmfunctions.runOpenMM(self, pdbPath, pdb, directory_dict["openmm_output"], directory_dict["intermediates"], self.params['forcefield'][0], self.params['implicit_solvent_model'][0])
        #else:
        #    for pdb in os.listdir(directory_dict["solvation_boxes"]):
        #        pdbPath = os.path.join(directory_dict["solvation_boxes"], pdb)
        #        self.log(f'Starting openMM with file: {pdbPath}')
        #        openmmfunctions.runOpenMM(self, pdbPath, pdb, directory_dict["md_raw_output"], directory_dict["intermediates"], self.params['forcefield'][0], self.params['explicit_solvent_model'][0])
        #    # Remove water molecules if implicit_solvent is False.
        #    for pdb in os.listdir(directory_dict["md_raw_output"]):
        #        pdbPath = os.path.join(directory_dict["md_raw_output"], pdb)
        #        openmmfunctions.removeWaterMolecules(pdbPath, pdb, directory_dict["openmm_output"])


        # export models
        for pdb in os.listdir(directory_dict["openmm_output"]):
            pdbPath = os.path.join(directory_dict["openmm_output"], pdb)
            pdbToExport = PDBFile(pdbPath)
            models_to_export.append(pdbToExport)

        self.output_models = models_to_export
        self.export_output_models()

        self.log('Completed OpenMM module run.')
        self.log('When you want to continue the haddock3 workflow after the OpenMM module, the next module should be topoaa, to rebuild the molecular topologies.')