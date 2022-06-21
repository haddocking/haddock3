"""OpenMM refinement module for HADDOCK3."""
from haddock.libs.libontology import Format, PDBFile
from haddock.modules import BaseHaddockModule
import os
import subprocess
from pathlib import Path
import haddock.modules.refinement.openmm.openmmfunctions as openmmfunctions

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")

class HaddockModule(BaseHaddockModule):
    """HADDOCK3 OpenMM module."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)
    
    def createDirectories(self):
        currentworkingdirectory = os.getcwd()
        openmmPdbfixer_output_Directory = os.path.join(currentworkingdirectory, 'openmm_pdbfixer_pdbs')
        modeller_solvationbox_pdbs_Directory = os.path.join(currentworkingdirectory, 'modeller_solvationbox_pdbs_Directory')
        openmm_intermediate_structures = os.path.join(currentworkingdirectory, 'openmm_intermediate_structures')
        openmm_md_raw_output_Directory = os.path.join(currentworkingdirectory, 'openmm_md_raw_output_Directory')
        openmm_output = os.path.join(currentworkingdirectory, 'openmm_output')
        self.log(f'Making directory for the pdbs builded with openmm pdbfixer: {openmmPdbfixer_output_Directory}')
        os.mkdir(openmmPdbfixer_output_Directory)
        self.log(f'Making directory for the pdbs with OpenMM modeller builded solvation boxes: {modeller_solvationbox_pdbs_Directory}')
        os.mkdir(modeller_solvationbox_pdbs_Directory)
        self.log(f'Making directory for saving the intermediate structures if the parameter save_intermediate_simulation_structures is specified: {openmm_intermediate_structures}')
        os.mkdir(openmm_intermediate_structures)
        self.log(f'Making directory for the openmm raw output containing water molecules: {openmm_md_raw_output_Directory}')
        os.mkdir(openmm_md_raw_output_Directory)
        self.log(f'Making directory for the openmm output: {openmm_output}')
        os.mkdir(openmm_output)
        return [openmmPdbfixer_output_Directory, modeller_solvationbox_pdbs_Directory, openmm_intermediate_structures, openmm_md_raw_output_Directory, openmm_output]
    
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
        models_to_export = []

        openmmPdbfixer_output_Directory, modeller_solvationbox_pdbs_Directory, openmm_intermediate_structures, openmm_md_raw_output_Directory, openmm_output = self.createDirectories()

        for pdb in previous_models:
            if(pdb.file_type is Format.PDB):
                pdb_filePath = f'{Path(pdb.path) / pdb.file_name}'
                pdb_filepath_openmm_pdbfixer_directory = os.path.join(openmmPdbfixer_output_Directory, pdb.file_name)                
                self.log('Fixing pdb with openmm pdbfixer.')
                openmmfunctions.OpenmmPdbfixer(self, pdb_filePath, pdb.file_name, openmmPdbfixer_output_Directory)
                if(self.params['implicit_solvent'] == False):
                    self.log(f'Building solvation box for file: {pdb.file_name}')
                    contains_xray_crystallography_cell_data = openmmfunctions.Does_pdb_contain_xray_crystallography_cell_data(pdb_filePath)
                    openmmfunctions.createSolvationBox(self, pdb_filepath_openmm_pdbfixer_directory, pdb.file_name, modeller_solvationbox_pdbs_Directory, self.params['forcefield'][0], self.params['explicit_solvent_model'][0], contains_xray_crystallography_cell_data)

        self.log(f'Starting openMM simulations.')
        if(self.params['implicit_solvent']):
            for pdb in os.listdir(openmmPdbfixer_output_Directory):
                pdbPath = os.path.join(openmmPdbfixer_output_Directory, pdb)
                self.log(f'Starting openMM with file: {pdbPath}')
                openmmfunctions.runOpenMM(self, pdbPath, pdb, openmm_output, openmm_intermediate_structures, self.params['forcefield'][0], self.params['implicit_solvent_model'][0])
        else:
            for pdb in os.listdir(modeller_solvationbox_pdbs_Directory):
                pdbPath = os.path.join(modeller_solvationbox_pdbs_Directory, pdb)
                self.log(f'Starting openMM with file: {pdbPath}')
                openmmfunctions.runOpenMM(self, pdbPath, pdb, openmm_md_raw_output_Directory, openmm_intermediate_structures, self.params['forcefield'][0], self.params['explicit_solvent_model'][0])
        
        # Remove water molecules if implicit_solvent is False.
        if(self.params['implicit_solvent'] == False):
            for pdb in os.listdir(openmm_md_raw_output_Directory):
                pdbPath = os.path.join(openmm_md_raw_output_Directory, pdb)
                openmmfunctions.removeWaterMolecules(pdbPath, pdb, openmm_output)

        for pdb in os.listdir(openmm_output):
            pdbPath = os.path.join(openmm_output, pdb)
            pdbToExport = PDBFile(pdbPath)
            models_to_export.append(pdbToExport)

        self.output_models = models_to_export
        self.export_output_models()

        self.log('Completed OpenMM module run.')
        self.log('When you want to continue the haddock3 workflow after the OpenMM module, the next module should be topoaa, to rebuild the molecular topologies.')