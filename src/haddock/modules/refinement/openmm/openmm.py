"""Refinement module using OpenMM.

The potential of OpenMM can be exploited to perform potentially different
 tasks, such as:
- refine the models obtained from a thorough docking run;
- refine the models in the middle of a docking run. For example, it can be used
 to refine the models coming from a `rigidbody` module before `flexref` is
 executed;
- refine some molecules prior to their use in a thorough docking run.

To get a list of all possible parameters, run::
     haddock3-cfg -m openmm
     
This module will refine all models coming from the previous workflow step and
 send them to the next step in the workflow. If you want to use other modules
such as `flexref` or `emref` after the OpenMM module, you need to recreate the
topologies by simply adding a `[topoaa]` step in the workflow.
See examples in `examples/thirdparty/openmm` folder.
"""
import os
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.app import (
    PME,
    ForceField,
    Modeller,
    PDBReporter,
    PDBxFile,
    Simulation,
    statedatareporter,
    )
from openmm.app.pdbfile import PDBFile as openmmpdbfile
from openmm.unit import kelvin, nanometer, nanometers, picosecond, picoseconds, atmosphere, molar  # noqa: E501
from pdbfixer import PDBFixer
from pdbtools import pdb_delhetatm

from haddock import log
from haddock.libs.libontology import PDBFile


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
        self.output_filename = self.model.file_name.replace(".pdb", "_omm.pdb")

    def import_constraints(self):
        """Import right openmm constraints."""
        if self.params["constraints"] == "HBonds":
            from openmm.app import HBonds as constraints
        elif self.params["constraints"] == "AllBonds":
            from openmm.app import AllBonds as constraints
        elif self.params["constraints"] == "HAngles":
            from openmm.app import HAngles as constraints
        return constraints

    def get_pdb_filepath(self, folder=None):
        """
        Get correct path to pdb file.
        
        Parameters
        ----------
        folder : str
            Path to folder
        """
        if folder:
            pdb_filepath = os.path.join(self.directory_dict["pdbfixer"],
                                        self.model.file_name)
        else:
            if isinstance(self.model, PDBFile):
                pdb_filepath = f'{self.model.rel_path}'
            else:
                pdb_filepath = str(Path(self.model.path, self.model.file_name))
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
        """
        Create solvation box for an explicit solvent simulation.
        
        Parameters
        ----------
        solvent_model : str
            One of the OpenMM's solvent models.
        """
        log.info(f'Building solvation box for file: {self.model.file_name}')
        pdb_filepath = self.get_pdb_filepath(self.directory_dict["pdbfixer"])
        forcefield = self.params["forcefield"]
        try:
            if ('.cif' in pdb_filepath):
                openmmpdb = PDBxFile(pdb_filepath)
            else:
                openmmpdb = openmmpdbfile(pdb_filepath)
            # creating Modeller
            modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
            usedForcefield = ForceField(forcefield, solvent_model)
            modeller.deleteWater()
            modeller.addHydrogens(usedForcefield, pH=7.0)
            # check the existence of xray cell data
            padding = self.params['padding'] * nanometers
            ions_conc = self.params['ion_concentration'] * molar
            # with neutralize=False it doesn't neutralize
            log.info(f"Solvating and adding ions to the system:"
                     f"padding {padding}, ions conc. : {ions_conc}")
            modeller.addSolvent(usedForcefield,
                                padding=padding,
                                neutralize=True,
                                ionicStrength=ions_conc)
            # Add required extra particles for forcefield,
            # e.g. Drude particles.
            if self.params['add_extra_particles_for_forcefield']:
                log.info("adding extra particles")
                modeller.addExtraParticles(forcefield)
            # write solvation box
            box_path = os.path.join(self.directory_dict["solvation_boxes"],
                                    self.model.file_name
                                    )
            openmmpdbfile.writeFile(modeller.topology,
                                    modeller.positions,
                                    open(box_path, 'w')
                                    )

        except Exception as e:
            log.info(f"BUILDING SOLVATION BOX ERROR: {e}")

    def remove_water_and_ions(self):
        """Remove water from the output of an explicit solvent run."""
        log.info(f"Removing water and ions from {self.model.file_name}")
        input_filepath = os.path.join(self.directory_dict["md_raw_output"],
                                      self.output_filename
                                      )
        output_filepath = os.path.join(self.directory_dict["openmm_output"],
                                       self.output_filename
                                       )
        with open(input_filepath, 'r') as fh:
            pdbYieldedLines = pdb_delhetatm.run(fh)
            with open(output_filepath, 'w') as writefh:
                writefh.writelines(pdbYieldedLines)
        
    def run_openmm(self, inputPDBfile, output_directory, solvent_model):
        """
        Run openmm simulation of the model pdb.

        Parameters
        ----------
        inputPDBfile : str
            Path to the PDB file.

        output_folder : str
            Path to the output directory.

        solvent_model : str
            One of the OpenMM's solvent models.
        """
        pdb = openmmpdbfile(inputPDBfile)
        forcefield = ForceField(self.params['forcefield'], solvent_model)
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
        simulation.reporters.append(
            statedatareporter.StateDataReporter(
                f"simulation_stats/observables_{self.identificator}.dat",
                50,
                step=True,
                totalEnergy=True,
                potentialEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                )
            )

        platform = simulation.context.getPlatform()
        log.info(f"Running on platform {platform.getName()}")
        log.info(f"Estimated platform speed: {platform.getSpeed()}")
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        log.info(f"Running NVT equilibration for {self.model.file_name}")

        # warming up the system with a NVT equilibration
        eq_steps = self.params['equilibration_timesteps']
        log.info(f"Warming up the system to {self.params['temperature_kelvin']}"
                 f" K in {eq_steps} steps")
        nvt_sims = 10
        delta_temp = self.params['temperature_kelvin'] / nvt_sims
        for n in range(nvt_sims):
            temperature = (delta_temp + (n * delta_temp)) * kelvin
            integrator.setTemperature(temperature)
            simulation.step(int(eq_steps / nvt_sims))
        
        # NPT simulation:
        system.addForce(
            MonteCarloBarostat(
                1 * atmosphere,
                self.params['temperature_kelvin'] * kelvin
                )
            )
        simulation.context.reinitialize(True)
        output_filepath = os.path.join(output_directory, self.output_filename)
        log.info(f"output file: {output_filepath}")
        simulation.reporters.append(PDBReporter(output_filepath,
                                    self.params['simulation_timesteps'])
                                    )
        # save intermediates
        if self.params['save_intermediate']:
            int_filepath = os.path.join(self.directory_dict["intermediates"],
                                        self.model.file_name
                                        )
            # get number of steps between each intermediate saved conf.
            steps_intermediate = self.params['simulation_timesteps'] // \
                self.params['save_intermediate']
            simulation.reporters.append(
                PDBReporter(int_filepath,
                            steps_intermediate)
                )
        # running real simulation
        log.info(f"Running simulation for {self.model.file_name}")
        simulation.step(self.params['simulation_timesteps'])
        log.info(f'OpenMM simulation successful for: {output_filepath}')

    def run(self):
        """Run openmm simulation."""
        # using pdbfixer
        self.openmm_pdbfixer()
        # explicit solvent simulation
        if not self.params['implicit_solvent']:
            # define solvent model and create solvation box
            solvent_model = self.params['explicit_solvent_model']
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
            solvent_model = self.params['implicit_solvent_model']
            solvent = "implicit solvent"

        log.info(f'starting {solvent} openMM simulation with file: {pdbPath}')
        self.run_openmm(pdbPath, output_folder, solvent_model)

        # Remove water molecules if implicit_solvent is False.
        if not self.params['implicit_solvent']:
            self.remove_water_and_ions()
