from haddock.libs.libontology import PDBFile
from openmm import *
from openmm.app.pdbfile import PDBFile as openmmpdbfile
from openmm.app import *
from openmm.unit import *
import os
#from openmm import Vec3, LangevinMiddleIntegrator
#from openmm.app import (Modeller, PDBxFile, ForceField, PDBReporter, Simulation, PME)
#from openmm.unit import nanometers, nanometer, picosecond, picoseconds, kelvin
from pdbfixer import PDBFixer
from pdbtools import pdb_delresname
from pathlib import Path
from haddock import log

def moveAtoms(atomPositionList, nanometerBoxSize, addition=True):
    shiftedAtomPostions = []
    if addition:
        sign = 1.0
    else:
        sign = -1.0
    
    shift = sign * nanometerBoxSize / 2 * 10
    for Vector3 in atomPositionList:
        newAtomVec3PositionManual = Vec3(Vector3.x * 10 + shift, Vector3.y * 10 + shift, Vector3.z * 10 + shift)
        shiftedAtomPostions.append(newAtomVec3PositionManual)

    return shiftedAtomPostions


class OPENMM:
    """CAPRI class."""

    def __init__(
            self,
            identificator,
            model,
            path,
            directory_dict,
            params,
            ):
        """
        Initialize the class.

        Parameters
        ----------
        identificator : str
            The identificator of the object.
        model : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            The model to be evaluated.
        path : Path
            Reference that defines where output should be saved.
        directory_dict : dict
            Dictionary of the directories handled by the openmm module.
        params : dict
            The openmm parameters.
        """
        self.identificator = identificator
        self.model = model
        self.path = path
        self.directory_dict = directory_dict
        self.params = params
        # other parameters
        self.output = Path("bu")


    def contains_xray_cell_data(self):
        """Determines if the model pdb contains xray cell data."""
        pdb_filepath = f'{Path(self.model.path) / self.model.file_name}'
        openmmpdb = openmmpdbfile(pdb_filepath)
        modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
        if(modeller.topology.getUnitCellDimensions() is None):
            xray_cell_data = False
        else:
            xray_cell_data =  True
        return xray_cell_data


    def OpenmmPdbfixer(self):
        """Calls Pdbfixer on the model pdb."""
        pdb_filepath = f'{Path(self.model.path) / self.model.file_name}'
        log.info(f'Fixing pdb: {self.model.file_name}')
        # creating PDBFixer
        fixer = PDBFixer(filename=pdb_filepath)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        # calling PDBFixer
        output_filepath = os.path.join(self.directory_dict["pdbfixer"], self.model.file_name)
        openmmpdbfile.writeFile(fixer.topology, fixer.positions, open(output_filepath, 'w'))

    def create_solvation_box(self, solvent_model):
        """Creates solvation box for an explicit solvent simulation."""
        pdb_filepath = os.path.join(self.directory_dict["pdbfixer"], self.model.file_name)
        forcefield = self.params["forcefield"][0]
        try:
            if('.cif' in pdb_filepath):
                openmmpdb = PDBxFile(pdb_filepath)
            else:
                openmmpdb = openmmpdbfile(pdb_filepath)
            # creating Modeller
            modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
            usedForcefield = ForceField(forcefield, solvent_model)
            modeller.deleteWater()
            modeller.addHydrogens(usedForcefield, pH=7.0)
            # check the existence of xray cell data
            if self.xray_cell_data:
                modeller.addSolvent(usedForcefield, padding=self.params['xray_cell_padding']*nanometers) # with neutralize=False it doesn't neutralize
            else:
                boxsizeIn_nm = self.params['solvent_boxsize_nm']
                modeller.addSolvent(usedForcefield, boxSize=Vec3(boxsizeIn_nm, boxsizeIn_nm, boxsizeIn_nm)*nanometers, neutralize=False)

            # Add required extra particles for forcefield, e.g. Drude particles.
            if(self.params['add_extra_particles_for_forcefield']):
                log.info(f"adding extra particles")
                modeller.addExtraParticles(forcefield)
            # check centering atoms
            if(self.params['move_atoms_to_solvationbox_center']):
                log.info(f"centering atoms")
                AtomPositions = moveAtoms(modeller.positions, 10)
            else:
                log.info(f"centering of atoms is not performed")
                AtomPositions = modeller.positions
            # write solvation box
            solvation_box_path = os.path.join(self.directory_dict["solvation_boxes"], self.model.file_name)
            openmmpdbfile.writeFile(modeller.topology, AtomPositions, open(solvation_box_path, 'w'))

        except Exception as e:
            log.info(f"BUILDING SOLVATION BOX ERROR: {e}")

    def remove_water_molecules(self):
        """Removes water molecules from the output of an explicit solvent run."""
        input_filepath = os.path.join(
                            self.directory_dict["md_raw_output"],
                            self.model.file_name)
        output_filepath = os.path.join(
                            self.directory_dict["openmm_output"],
                            self.model.file_name)
        with open(input_filepath, 'r') as fh:
            pdbYieldedLines = pdb_delresname.run(fh, 'HOH')
            with open(output_filepath, 'w') as writefh:
                writefh.writelines(pdbYieldedLines)

    def runOpenMM(self, inputPDBfile, output_directory, solvent_model):
        """Runs openmm simulation of the model pdb."""
        pdb = openmmpdbfile(inputPDBfile)
        forcefield = ForceField(self.params['forcefield'][0], solvent_model)
        system = forcefield.createSystem(pdb.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer,
            constraints=globals()[self.params['constraints']],
            removeCMMotion=self.params['remove_center_of_mass_motion'],
            rigidWater=self.params['rigid_water']) #ERROR when system_constraints = 'None'.
        integrator = LangevinMiddleIntegrator(
            self.params['temperature_kelvin']*kelvin,
            1/picosecond,
            self.params['timestep_ps']*picoseconds
            )
        # integrator seed
        if self.params["seed"]:
            log.info(f"setting non-random seed = {self.params['seed']}") #set_seed = 42
            integrator.setRandomNumberSeed(self.params["seed"])
        
        get_seed = integrator.getRandomNumberSeed()
        log.info(f"simulation seed {get_seed}")
        if self.params["seed"] != 0:
            if get_seed != self.params["seed"]:
                log.warning(f"integrator seed ({get_seed}) does not match input seed ({self.params['seed']})")
        # simulation
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.step(self.params['equilibration_timesteps'])
        output_filepath = os.path.join(output_directory, self.model.file_name)
        simulation.reporters.append(PDBReporter(output_filepath, self.params['simulation_timesteps']))
        # save intermediates
        if(self.params['save_intermediate']):
            int_filepath = os.path.join(self.directory_dict["intermediates"], self.model.file_name)
            simulation.reporters.append(
                PDBReporter(int_filepath,
                            self.params['steps_intermediate'])
                            )
        # simulation.reporters.append(StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True, temperature=True)) # Report system state of the simulation.
        simulation.step(self.params['simulation_timesteps'])
        log.info(f'OpenMM simulation succesful for: {output_filepath}')

    def run(self):
        """Run openmm simulation."""
        # using pdbfixer
        self.OpenmmPdbfixer()
        # explicit solvent simulation
        if(self.params['implicit_solvent'] == False):
            # solvation box
            log.info(f'Building solvation box for file: {self.model.file_name}')
            # check if pdb contains xray cell data
            self.xray_cell_data = self.contains_xray_cell_data()
            solvent_model = self.params['explicit_solvent_model'][0]
            self.create_solvation_box(solvent_model)
            # output files
            pdbPath = os.path.join(self.directory_dict["solvation_boxes"], self.model.file_name)
            output_folder = self.directory_dict["md_raw_output"]
        # implicit solvent simulation
        else: 
            pdbPath = os.path.join(self.directory_dict["pdbfixer"], self.model.file_name)
            output_folder = self.directory_dict["openmm_output"]
            solvent_model = self.params['implicit_solvent_model'][0]
        
        log.info(f'starting openMM simulation with file: {pdbPath}')
        self.runOpenMM(pdbPath, output_folder, solvent_model)
        
        # Remove water molecules if implicit_solvent is False.
        if self.params['implicit_solvent'] == False:
            self.remove_water_molecules()
        