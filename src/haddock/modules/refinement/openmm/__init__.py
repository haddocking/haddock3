"""OpenMM refinement module for HADDOCK3."""
import subprocess
from pathlib import Path

from haddock.libs.libontology import Format, PDBFile
from haddock.modules import BaseHaddockModule
from haddock.modules.refinement.openmm import openmmfunctions


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 OpenMM module."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

    def create_directories(self):
        """Create directories required for OpenMM."""
        cwd = Path.cwd()

        self.openmmPdbfixer_output_Directory = Path(cwd, 'openmm_pdbfixer_pdbs')

        self.modeller_solvationbox_pdbs_Directory = \
            Path(cwd, 'self.modeller_solvationbox_pdbs_Directory')

        self.openmm_intermediate_structures = \
            Path(cwd, 'self.openmm_intermediate_structures')

        self.openmm_md_raw_output_Directory = \
            Path(cwd, 'self.openmm_md_raw_output_Directory')

        self.openmm_output = Path(cwd, 'self.openmm_output')

        self.log(
            'Making directory for the pdbs builded with openmm pdbfixer: '
            f'{self.openmmPdbfixer_output_Directory}'
            )
        self.openmmPdbfixer_output_Directory.mkdir()

        self.log(
            'Making directory for the pdbs with OpenMM modeller builded '
            f'solvation boxes: {self.modeller_solvationbox_pdbs_Directory}')
        self.modeller_solvationbox_pdbs_Directory.mkdir()

        self.log(
            'Making directory for saving the intermediate structures if the '
            'parameter save_intermediate_simulation_structures is specified: '
            f'{self.openmm_intermediate_structures}'
            )
        self.openmm_intermediate_structures.mkdir()

        self.log(
            'Making directory for the openmm raw output containing water '
            f'molecules: {self.openmm_md_raw_output_Directory}'
            )
        self.openmm_md_raw_output_Directory.mkdir()

        self.log(
            'Making directory for the openmm output: '
            f'{self.openmm_output}'
            )
        self.openmm_output.mkdir()

        return

    @classmethod
    def confirm_installation(cls):
        """Confirm OpenMM is properly installed."""
        def runSubprocess(command_to_run):
            output = subprocess.run(
                [command_to_run],
                shell=True,
                capture_output=True,
                encoding='utf-8',
                )
            return output.stdout.strip()

        output_openmm = runSubprocess("conda list openmm --json")
        output_pdbfixed = runSubprocess("conda list pdbfixer --json")

        if(output_openmm == '[]'):
            raise Exception('OpenMM is not installed in conda.')

        if(output_pdbfixed == '[]'):
            raise Exception('OpenMM pdbfixer is not installed in conda.')
        return

    def _run(self):
        """Execute module."""
        previous_models = self.previous_io.retrieve_models()[0]
        models_to_export = []

        self.create_directories()

        for pdb in previous_models:
            if pdb.file_type is Format.PDB:

                pdb_filePath = Path(pdb.path, pdb.file_name)

                self.pdb_filepath_openmm_pdbfixer_directory = \
                    Path(self.openmmPdbfixer_output_Directory, pdb.file_name)

                self.log('Fixing PDB with OpenMM PDBFixer.')

                openmmfunctions.OpenmmPdbfixer(
                    self,
                    pdb_filePath,
                    pdb.file_name,
                    self.openmmPdbfixer_output_Directory,
                    )

                if not self.params['implicit_solvent']:
                    self.log(
                        'Building solvation box for file: '
                        f'{pdb.file_name}'
                        )

                    contains_xray_crystallography_cell_data = \
                        openmmfunctions.Does_pdb_contain_xray_crystallography_cell_data(pdb_filePath)  # noqa: E501

                    openmmfunctions.createSolvationBox(
                        self,
                        self.pdb_filepath_openmm_pdbfixer_directory,
                        pdb.file_name,
                        self.modeller_solvationbox_pdbs_Directory,
                        self.params['forcefield'][0],
                        self.params['explicit_solvent_model'][0],
                        contains_xray_crystallography_cell_data,
                        )

        self.log('Starting openMM simulations.')
        if self.params['implicit_solvent']:

            for pdb_path in self.openmmPdbfixer_output_Directory.iterdir():
                self.log(f'Starting openMM with file: {str(pdb_path)}')
                openmmfunctions.runOpenMM(
                    self,
                    pdb_path,
                    pdb_path.name,
                    self.openmm_output,
                    self.openmm_intermediate_structures,
                    self.params['forcefield'][0],
                    self.params['implicit_solvent_model'][0],
                    )

        else:
            for pdb_path in self.modeller_solvationbox_pdbs_Directory.iterdir():
                self.log(f'Starting openMM with file: {str(pdb_path)}')
                openmmfunctions.runOpenMM(
                    self,
                    pdb_path,
                    pdb_path.name,
                    self.openmm_md_raw_output_Directory,
                    self.openmm_intermediate_structures,
                    self.params['forcefield'][0],
                    self.params['explicit_solvent_model'][0],
                    )

        # Remove water molecules if implicit_solvent is False.
        if not self.params['implicit_solvent']:

            for pdb_path in self.openmm_md_raw_output_Directory.iterdir():
                openmmfunctions.removeWaterMolecules(
                    pdb_path,
                    pdb_path.name,
                    self.openmm_output,
                    )

        for pdb_path in self.openmm_output.iterdir():
            pdb_to_export = PDBFile(pdb_path)
            models_to_export.append(pdb_to_export)

        self.output_models = models_to_export
        self.export_output_models()

        self.log('Completed OpenMM module run.')
        self.log(
            'When you want to continue the haddock3 workflow after the OpenMM '
            'module, the next module should be topoaa, '
            'to rebuild the molecular topologies.'
            )
