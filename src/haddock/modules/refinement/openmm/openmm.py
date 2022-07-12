"""OpenMM module."""
import os
from pathlib import Path

from openmm import LangevinMiddleIntegrator, Vec3
from openmm.app import (
    PME,
    ForceField,
    Modeller,
    PDBReporter,
    PDBxFile,
    Simulation,
    )
from openmm.app.pdbfile import PDBFile as openmmpdbfile
from openmm.unit import kelvin, nanometer, nanometers, picosecond, picoseconds
from pdbfixer import PDBFixer
from pdbtools import pdb_delresname

from haddock import log
from haddock.libs.libontology import PDBFile


def moveAtoms(atomPositionList, nanometerBoxSize, addition=True):
    """Shift atoms positions."""
    shiftedAtomPostions = []
    if addition:
        sign = 1.0
    else:
        sign = -1.0
    
    shift = sign * nanometerBoxSize / 2 * 10
    for Vector3 in atomPositionList:
        newAtomVec3PositionManual = Vec3(Vector3.x * 10 + shift,
                                         Vector3.y * 10 + shift,
                                         Vector3.z * 10 + shift
                                         )
        shiftedAtomPostions.append(newAtomVec3PositionManual)

    return shiftedAtomPostions


class OPENMM:
    """OPENMM class."""

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
        self.output = Path("output_openmm.log")
        self.constraints = self.import_constraints()

    def import_constraints(self):
        """Import right openmm constraints."""
        if self.params["constraints"] == "HBonds":
            from openmm.app import HBonds as constraints
        elif self.params["constraints"] == "AllBonds":
            from openmm.app import AllBonds as constraints
        elif self.params["constraints"] == "HAngles":
            from openmm.app import HAngles as constraints
        return constraints

    def contains_xray_cell_data(self):
        """Determine if the model pdb contains xray cell data."""
        pdb_filepath = self.get_pdb_filepath()
        openmmpdb = openmmpdbfile(pdb_filepath)
        modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
        if(modeller.topology.getUnitCellDimensions() is None):
            xray_cell_data = False
        else:
            xray_cell_data = True
        return xray_cell_data
    
    def get_pdb_filepath(self, folder=None):
        """Get correct path to pdb file."""
        if folder:
            pdb_filepath = os.path.join(self.directory_dict["pdbfixer"],
                                        self.model.file_name)
        else:
            if isinstance(self.model, PDBFile):
                pdb_filepath = f'{self.model.rel_path}'
            else:
                pdb_filepath = f'{Path(self.model.path) / self.model.file_name}'
        log.debug(f"pdb_filepath is {pdb_filepath}")
        return pdb_filepath

    def openmm_pdbfixer(self):
        """Call Pdbfixer on the model pdb."""
        pdb_filepath = self.get_pdb_filepath()
        log.info(f'Fixing pdb: {self.model.file_name} (path {pdb_filepath})')
        # creating PDBFixer
        fixer = PDBFixer(filename=pdb_filepath)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        # calling PDBFixer
        output_filepath = os.path.join(self.directory_dict["pdbfixer"],
                                       self.model.file_name
                                       )
        openmmpdbfile.writeFile(fixer.topology,
                                fixer.positions,
                                open(output_filepath, 'w')
                                )

    def create_solvation_box(self, solvent_model):
        """Create solvation box for an explicit solvent simulation."""
        pdb_filepath = self.get_pdb_filepath(self.directory_dict["pdbfixer"])
        forcefield = self.params["forcefield_fname"]
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
                xray_padding = self.params['xray_cell_padding'] * nanometers
                # with neutralize=False it doesn't neutralize
                modeller.addSolvent(usedForcefield,
                                    padding=xray_padding)
            else:
                boxsizeIn_nm = self.params['solvent_boxsize_nm']
                box_size = Vec3(boxsizeIn_nm,
                                boxsizeIn_nm,
                                boxsizeIn_nm) * nanometers
                modeller.addSolvent(usedForcefield,
                                    boxSize=box_size,
                                    neutralize=False
                                    )

            # Add required extra particles for forcefield,
            # e.g. Drude particles.
            if(self.params['add_extra_particles_for_forcefield']):
                log.info("adding extra particles")
                modeller.addExtraParticles(forcefield)
            # check centering atoms
            if(self.params['move_atoms_to_solvationbox_center']):
                log.info("centering atoms")
                AtomPositions = moveAtoms(modeller.positions, 10)
            else:
                log.info("centering of atoms is not performed")
                AtomPositions = modeller.positions
            # write solvation box
            box_path = os.path.join(self.directory_dict["solvation_boxes"],
                                    self.model.file_name
                                    )
            openmmpdbfile.writeFile(modeller.topology,
                                    AtomPositions,
                                    open(box_path, 'w')
                                    )

        except Exception as e:
            log.info(f"BUILDING SOLVATION BOX ERROR: {e}")

    def remove_water_molecules(self):
        """Remove water from the output of an explicit solvent run."""
        input_filepath = os.path.join(self.directory_dict["md_raw_output"],
                                      self.model.file_name
                                      )
        output_filepath = os.path.join(self.directory_dict["openmm_output"],
                                       self.model.file_name
                                       )
        with open(input_filepath, 'r') as fh:
            pdbYieldedLines = pdb_delresname.run(fh, 'HOH')
            with open(output_filepath, 'w') as writefh:
                writefh.writelines(pdbYieldedLines)

    def runOpenMM(self, inputPDBfile, output_directory, solvent_model):
        """Run openmm simulation of the model pdb."""
        pdb = openmmpdbfile(inputPDBfile)
        forcefield = ForceField(self.params['forcefield_fname'], solvent_model)
        # system setup
        # should give an ERROR when system_constraints = 'None'.
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1 * nanometer,
            constraints=self.constraints,
            removeCMMotion=self.params['remove_center_of_mass_motion'],
            rigidWater=self.params['rigid_water']
            )
        # integrator definition. Few freedom here
        integrator = LangevinMiddleIntegrator(
            self.params['temperature_kelvin'] * kelvin,
            1 / picosecond,
            self.params['timestep_ps'] * picoseconds
            )
        
        # setting integrator seed
        seed = self.params["seed"]
        if seed != 0:
            log.info(f"setting non-random seed = {self.params['seed']}")
            integrator.setRandomNumberSeed(self.params["seed"])
        # retrieving integrator seed
        get_seed = integrator.getRandomNumberSeed()
        log.info(f"simulation seed {get_seed} for {self.model.file_name}")

        # simulation
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        log.info(f"Running equilibration for {self.model.file_name}")
        simulation.step(self.params['equilibration_timesteps'])
        output_filepath = os.path.join(output_directory, self.model.file_name)
        simulation.reporters.append(PDBReporter(output_filepath,
                                    self.params['simulation_timesteps'])
                                    )
        # save intermediates
        if self.params['save_intermediate']:
            int_filepath = os.path.join(self.directory_dict["intermediates"],
                                        self.model.file_name
                                        )
            simulation.reporters.append(
                PDBReporter(int_filepath,
                            self.params['steps_intermediate'])
                )
        # simulation.reporters.append(StateDataReporter(sys.stdout,
        # 100,
        # step=True,
        # potentialEnergy=True,
        # temperature=True)
        # )
        # # Report system state of the simulation.
        simulation.step(self.params['simulation_timesteps'])
        log.info(f'OpenMM simulation succesful for: {output_filepath}')

    def run(self):
        """Run openmm simulation."""
        # using pdbfixer
        self.openmm_pdbfixer()
        # explicit solvent simulation
        if not self.params['implicit_solvent']:
            # solvation box
            log.info(f'Building solvation box for file: {self.model.file_name}')
            # check if pdb contains xray cell data
            self.xray_cell_data = self.contains_xray_cell_data()
            solvent_model = self.params['explicit_solvent_model_fname']
            self.create_solvation_box(solvent_model)
            # output files and string
            pdbPath = os.path.join(self.directory_dict["solvation_boxes"],
                                   self.model.file_name
                                   )
            output_folder = self.directory_dict["md_raw_output"]
            solvent = "explicit solvent"
        # implicit solvent simulation
        else:
            pdbPath = os.path.join(self.directory_dict["pdbfixer"],
                                   self.model.file_name
                                   )
            output_folder = self.directory_dict["openmm_output"]
            solvent_model = self.params['implicit_solvent_model_fname']
            solvent = "implicit solvent"
        
        log.info(f'starting {solvent} openMM simulation with file: {pdbPath}')
        self.runOpenMM(pdbPath, output_folder, solvent_model)
        
        # Remove water molecules if implicit_solvent is False.
        if not self.params['implicit_solvent']:
            self.remove_water_molecules()
