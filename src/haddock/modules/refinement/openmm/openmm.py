"""
OpenMM Molecular Dynamics refinement module for HADDOCK3.
"""

# Standard libarires importation
import os

from contextlib import suppress
from pathlib import Path

# Import OpenMM and pdbfixer third-party libraries
# >conda activate haddock3
# >conda install -c conda-forge libstdcxx-ng openmm pdbfixer

# allow general testing when OpenMM is not installed
#with suppress(ImportError):
from openmm import (
    CustomCentroidBondForce,
    CustomExternalForce,
    LangevinMiddleIntegrator,
    MonteCarloBarostat,
    System,
    )
from openmm.app import (
    PME,
    NoCutoff,
    AllBonds,
    ForceField,
    HAngles,
    HBonds,
    Modeller,
    PDBxFile,
    Simulation,
    statedatareporter,
    )
from openmm.app.pdbfile import PDBFile as openmmpdbfile
from openmm.unit import (
    MOLAR_GAS_CONSTANT_R,
    angstroms,
    atmosphere,
    femtoseconds,
    kelvin,
    kilocalorie_per_mole,
    molar,
    nanometer,
    nanometers,
    picosecond,
    picoseconds,
    )
from pdbfixer import PDBFixer

# Haddock libraries
from pdbtools import pdb_delhetatm, pdb_mkensemble
from haddock import log
from haddock.core.exceptions import ModuleError
from haddock.core.typing import Optional, Union, ParamDict
from haddock.libs.libontology import PDBFile
from haddock.libs.libpdb import add_TER_on_chain_breaks


class OPENMM:
    """OPENMM class."""

    def __init__(
            self,
            identificator: int,
            model: PDBFile,
            path: Path,
            directory_dict: dict[str, str],
            params: ParamDict,
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
        directory_dict : dict[str, str]
            Dictionary of the directories handled by the openmm module.
            {d: d
                for d in [
                    "pdbfixer",
                    "solvation_boxes",
                    "intermediates",
                    "md_raw_output",
                    "openmm_output",
                    "simulation_stats"
                    ]
                }
        params : ParamDict
            The openmm modules parameters.
        """
        self.identificator = identificator
        self.model = model
        self.path = path
        self.directory_dict = directory_dict
        self.params = params
        self.log = log

        # other parameters
        self.output = Path("output_openmm.log")
        self.constraints = self.import_constraints()
        self.output_filename = self.model.file_name.replace(".pdb", "_omm.pdb")

    def import_constraints(self):  # type: ignore
        """Cast parameter string to proper openmm constraints."""
        if self.params["constraints"] == "HBonds":
            return HBonds
        elif self.params["constraints"] == "AllBonds":
            return AllBonds
        elif self.params["constraints"] == "HAngles":
            return HAngles
        return None

    def get_pdb_filepath(self, folder: Union[bool, str] = None) -> str:
        """Get correct path to pdb file.
        
        Parameters
        ----------
        folder : str
            Path to folder
        """
        if folder:
            pdb_filepath = os.path.join(
                self.directory_dict["pdbfixer"],
                self.model.file_name,
                )
        else:
            if isinstance(self.model, PDBFile):
                pdb_filepath = f"{self.model.rel_path}"
            else:
                pdb_filepath = str(Path(self.model.path, self.model.file_name))
        self.log.debug(f"pdb_filepath is {pdb_filepath}")
        return pdb_filepath
    
    def handle_chainbreaks(self, pdb_fpath: str) -> str:
        output_filepath = os.path.join(
            self.directory_dict["pdbfixer"],
            f"withTER_{self.model.file_name}",
            )
        add_TER_on_chain_breaks(pdb_fpath, output_filepath)
        return output_filepath

    def openmm_pdbfixer(self):
        """Call Pdbfixer on the model pdb."""
        input_pdb_filepath = self.get_pdb_filepath()
        pdb_filepath = self.handle_chainbreaks(input_pdb_filepath)
        self.log.info(
            f"Fixing pdb: {self.model.file_name} (path {pdb_filepath})"
            )
        # creating PDBFixer
        fixer = PDBFixer(filename=pdb_filepath)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(self.params["ph"])
        # calling PDBFixer
        output_filepath = os.path.join(
            self.directory_dict["pdbfixer"],
            self.model.file_name,
            )
        openmmpdbfile.writeFile(
            fixer.topology,
            fixer.positions,
            open(output_filepath, 'w'),
            keepIds=True,
            )

    def create_solvation_box(self, solvent_model: str) -> Optional[str]:
        """Create solvation box for an explicit solvent simulation.

        Parameters
        ----------
        solvent_model : str
            One of the OpenMM's solvent models.

        Returns
        -------
        box_path : Optional[str]
            Path to the solvated system

        Raise
        ----------
        Exception
            Not able to generate the solvation box.
        """
        self.log.info(f"Building solvation box for: {self.model.file_name}")
        pdb_filepath = self.get_pdb_filepath(self.directory_dict["pdbfixer"])
        forcefield = self.params["forcefield"]
        try:
            if ('.cif' in pdb_filepath):
                openmmpdb = PDBxFile(pdb_filepath)
            else:
                openmmpdb = openmmpdbfile(pdb_filepath)
            # Obtain initial input chains
            self.input_chains = [
                chain.id
                for chain in openmmpdb.topology.chains()
                ]
            # creating Modeller
            modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
            usedForcefield = ForceField(forcefield, solvent_model)
            modeller.deleteWater()
            modeller.addHydrogens(usedForcefield, pH=self.params["ph"])
            # check the existence of xray cell data
            padding = self.params["padding"] * nanometers
            ions_conc = self.params["ion_concentration"] * molar
            # with neutralize=False it doesn't neutralize
            log.info(
                f"Solvating and adding ions to the system:"
                f"padding {padding}, ions conc. : {ions_conc}"
                )
            modeller.addSolvent(
                usedForcefield,
                padding=padding,
                neutralize=True,
                ionicStrength=ions_conc,
                )
            # Add required extra particles for forcefield,
            # e.g. Drude particles.
            if self.params["add_extra_particles_for_forcefield"]:
                log.info("adding extra particles")
                modeller.addExtraParticles(forcefield)
            # write solvation box
            box_path = os.path.join(
                self.directory_dict["solvation_boxes"],
                self.model.file_name
                )
            openmmpdbfile.writeFile(
                modeller.topology,
                modeller.positions,
                open(box_path, "w"),
                keepIds=True,
                )

        except Exception as e:
            error_msg = (
                f"An error occured when building solvation box: {e}"
                )
            self.log.error(error_msg)
            raise ModuleError(f"[openMM] module error: {error_msg}")
        else:
            return box_path

    def equilibrate_solvation_box(
            self,
            pdb_filepath: str,
            solvent_model: str,
            ) -> Optional[str]:
        """Machinery for the equilibration of water in presence of the protein.

        Here, the idea is to:
            0. Check if something has to be done
            1. Initiate the system
            2. Construct restrain forces on protein atom positions
            3. Run simulation under restrain
            4. Save new coordinates file

        Parameters
        ----------
        pdb_filepath : str
            Path to a water boxed system
        solvent_model : str
            One of the OpenMM's solvent models.

        Returns
        -------
        eq_pdb_filepath: str
            Path to an equilibrated water boxed system, readdy of MD

        Raises
        ------
        Exception
            Unable to preform the equilibration.
        """
        # Check if spring constant > 0
        if not self.params["solv_equilibration"]:
            return pdb_filepath
        
        self.log.info(f"Equilibrating solvation box from {pdb_filepath}...")
        try:
            # 1. Initiate the system and simulation objects
            # Load structure coordinates
            system_coordinates = openmmpdbfile(pdb_filepath)
            # Initiate forcefiled object
            forcefield = ForceField(self.params["forcefield"], solvent_model)
            # system setup
            # should give an ERROR when system_constraints = 'None'.
            system = forcefield.createSystem(
                system_coordinates.topology,
                nonbondedMethod=PME,
                nonbondedCutoff=1 * nanometer,
                constraints=HAngles,
                removeCMMotion=self.params["remove_center_of_mass_motion"],
                rigidWater=self.params["rigid_water"],
                )
            # Initiate intergrator
            integrator = LangevinMiddleIntegrator(
                self.params["solv_eq_max_temperature_kelvin"] * kelvin,
                1 / picosecond,
                self.params["solv_eq_stepsize_fs"] * femtoseconds
                )
            # Set pseudo-random seed
            integrator.setRandomNumberSeed(self.params["iniseed"])

            # 2. Add constraints to the protein atoms
            self.log.info("Adding protein constraints")
            # Add position restraint force
            rest_force, spring_force_index = self._gen_restrain_force(
                system_coordinates.topology.atoms(),
                system_coordinates.positions,
                spring_constant=self.params["solv_eq_spring_constant"],
                )
            # Add the restrain force to the system
            system.addForce(rest_force)

            # Initiate simulation object
            simulation = Simulation(
                system_coordinates.topology,
                system,
                integrator,
                )

            # Initiate state reporter object
            reported_output_fname = (
                f'{self.directory_dict["simulation_stats"]}/'
                f"solv_eq_observables_{self.identificator}.dat"
                )
            statereporter = statedatareporter.StateDataReporter(
                reported_output_fname,
                10,
                step=True,
                totalEnergy=True,
                potentialEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                )
            simulation.reporters.append(statereporter)
            simulation.context.setPositions(system_coordinates.positions)

            # 3. Run simulation under constraints
            # Process energy minimization prior to any dynamic simulation
            self.log.info("Minimizing system energy")
            simulation.minimizeEnergy()
            # Set parameters
            eq_steps = self.params["solv_eq_timesteps"]
            nvt_sims = 7
            delta_steps = max(1, int(eq_steps / nvt_sims))
            max_temperature = self.params["solv_eq_max_temperature_kelvin"]
            max_temperature = max_temperature
            delta_temp = max_temperature / nvt_sims
            # Progressively warmup the system
            self.log.info("Warming system up...")
            integrator.setFriction(100 / picosecond)
            for n in range(1, nvt_sims + 1):
                # Compute new values
                qtemp = n * delta_temp
                # Set values
                integrator.setTemperature(qtemp * kelvin)
                # Run some md steps
                simulation.step(delta_steps)
            # Set final temperature and temperature coupling
            integrator.setTemperature(max_temperature * kelvin)
            integrator.setFriction(1 / picosecond)
            # Makes sure temperature of the system is reached
            self._stabilize_temperature(
                simulation,
                max_temperature,
                statereporter._dof,
                tolerance=10.0,
                steps=50
                )
            # Print log info
            self.log.info(
                "Running solvent equilibration phase "
                "under protein coordinate position restraints:"
                )
            self.log.info(
                f"For {self.params['solv_eq_timesteps']} timesteps at "
                f"{self.params['solv_eq_max_temperature_kelvin']} K "
                f"(stepsize={self.params['solv_eq_stepsize_fs']} fs)..."
                )
            # Do few simulation steps at max temperature
            # progressively releasing the position restraints on the protein
            slow_reduction_steps = 10
            red_delta_steps = int(eq_steps / slow_reduction_steps)
            for n in range(slow_reduction_steps - 1, -1, -1):
                const_stength = n / slow_reduction_steps
                new_spring_const = self.params["solv_eq_spring_constant"] * const_stength
                rest_force.setGlobalParameterDefaultValue(
                    spring_force_index,
                    new_spring_const,
                    )
                simulation.step(red_delta_steps)
            # Progressively freeze the system
            self.log.info("Cooling down the system...")
            for n in range(nvt_sims - 1, 0, -1):
                qtemp = n * delta_temp
                integrator.setTemperature(qtemp * kelvin)
                simulation.step(delta_steps)

            # 4. Save new coordinates file
            # Write finale frame structure
            eq_pdb_filepath = os.path.join(
                self.directory_dict["solvation_boxes"],
                f"solvent_eq_{self.model.file_name}",
                )
            # Retrieve last coordinates
            eq_state = simulation.context.getState(getPositions=True)
            eq_positions = eq_state.getPositions()
            openmmpdbfile.writeFile(
                simulation.topology,
                eq_positions,
                open(eq_pdb_filepath, "w"),
                keepIds=True,
                )
            self.log.info(f"Equilibrated solvated system: {eq_pdb_filepath} !")
            return eq_pdb_filepath

        except Exception:
            import traceback
            strerr = traceback.format_exc()
            self.log.debug(
                "Error durring solvent equilibration for file "
                f"{pdb_filepath}:\n{strerr}.{os.linesep}"
                "Continuing from unequilibrated solvent frame."
                )
            return pdb_filepath

    def _stabilize_temperature(
            self,
            simulation: Simulation,
            temperature: float,
            dof: int,
            tolerance: float = 5.0,
            steps: int = 50,
            ) -> int:
        """
        Make sure the simulated system reached the desired temperature.
        
        Parameters
        ----------
        simulation : py:class:`openmm.Simulation`
            An openmm system
        temperature : float
            The temperature hoped to be reached
        dof : int
            The degree of freedom obtained from the statereporter._dof
        tolerance : float
            The tolerance allowed for the temperature
        steps : int
            The number of steps to do before checking again that
            temperature was reached
        """
        additional_steps: int = 0
        # Makes sure temperature of the system is reached
        while not ((temperature - tolerance)
                   <= self._get_simulation_temperature(simulation, dof)
                   <= (temperature + tolerance)):
            # Do several simulation steps
            simulation.step(steps)
            additional_steps += steps
        return additional_steps

    @staticmethod
    def _gen_restrain_force(
            atoms: list,
            positions: list,
            subset: Union[bool, list] = None,
            spring_constant: float = 20.0,
            ) -> CustomExternalForce:
        """
        Generate CustomExternalForce aiming at restraining protein coordinates.

        Parameters
        ----------
        atoms : list
            List of the particules present in the topology of the system
        positions : list
            Initial position of the particules in the system
        subset : list
            List of atom indexes on which to apply the restrain force
        spring_constant : float
            Initial value of the position restrain spring constant.

        Return
        ----------
        rest_force: :py:class:`openmm.openmm.CustomExternalForce`
            Set CustomExternalForce object restraining protein position
        spring_force_index: int
            Index of the spring force constant
        """
        # Initiate Custom Force to restrain positions
        rest_force = CustomExternalForce(
            # could have been "k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
            # but because we are in periodic boundary conditions,
            # it is better to use "k*periodicdistance(x, y, z, x0, y0, z0)^2"
            # c.f. OpenMM documentation
            "k*periodicdistance(x, y, z, x0, y0, z0)^2"
            )
        # Define spring constant value
        spring_force_index = rest_force.addGlobalParameter(
            "k",
            spring_constant * kilocalorie_per_mole / angstroms ** 2
            )
        rest_force.addPerParticleParameter('x0')
        rest_force.addPerParticleParameter('y0')
        rest_force.addPerParticleParameter('z0')
        # Set of atom indexes to work on
        ind_subset = list(range(len(positions))) if not subset else subset
        # Loop over topology atoms
        for atom in atoms:
            # Filter out non-protein particles and hydrogens
            if atom.element.symbol == 'H':  # Hydrogen atom
                continue
            if atom.residue.name == 'HOH':  # Water molecule
                continue
            if len(atom.residue.name) != 3:  # Skip ions
                continue
            if atom.index not in ind_subset:  # Not in subset
                continue
            # Hold index of the restrained atom
            rest_force.addParticle(atom.index, positions[atom.index])
        return rest_force, spring_force_index

    def _gen_centroid_forces(self, system: System) -> None:
        """Add centroid forces between protein monomers to the system.

        FIXME: Not functionning due to issue of group definitions in
               the self._chain_centroid_force() function

        Parameters
        ----------
        system : py:class:`openmm.System`
            An openmm system
        """
        for i, chainid in enumerate(self.input_chains[:-1]):
            for chain2 in self.input_chains[i + 1:]:
                new_force = self._chain_centroid_force(system, chainid, chain2)
                system.addForce(new_force)

    @staticmethod
    def _get_chain_atoms(system: System, chainid: str) -> list[int]:
        """Retrun list of atom ids belonging to a chain.
        
        Parameters
        ----------
        system : py:class:`openmm.System`
            An openmm system
        chainid : str
            Chain id on which atom ids must be retrieved from
        """
        for chain in system.topology.chains():
            if chain.id == chainid:
                return [atom.id for atom in chain.atoms()]
        return []

    def _chain_centroid_force(
            self,
            system: System,
            chain1: str,
            chain2: str,
            ) -> CustomCentroidBondForce:
        """Set up a centroid force betweem two chains.

        FIXME: Not functionning ...

        Parameters
        ----------
        system : py:class:`openmm.System`
            An openmm system
        chain1 : str
            One of system chain (hopfully a monomer) to apply the force on
        chain2 : str
            Second chain (hopfully a monomer) to apply the force on
        """
        centroid_force = CustomCentroidBondForce(2, "0.1*k*distance(g1,g2)^2")
        centroid_force.addPerBondParameter("k")
        gc1 = centroid_force.addGroup(self._get_chain_atoms(system, chain1))
        gc2 = centroid_force.addGroup(self._get_chain_atoms(system, chain2))
        centroid_force.addBond([gc1, gc2], [1])
        return centroid_force

    @staticmethod
    def _get_simulation_temperature(simulation: Simulation, dof: int) -> float:
        """Computes/gather the current system temperature.

        NOTE: this function is similar to the lines of code found in
              the :py:class:`openmm.app.statedatareporter.StateDataReporter()`
        
        NOTE2: Would be better if in the StateustomCentroidBondForce class...

        Parameters
        ----------
        simulation : :py:class:`openmm.app.simulation.Simulation()`
            The current simulation
        dof : int
            Value of the StateDataReporter._dof

        Return
        ----------
        temperature_k : float
            The current simulation temperature (in kelvin)
        """
        integrator = simulation.context.getIntegrator()
        if hasattr(integrator, "computeSystemTemperature"):
            temperature = integrator.computeSystemTemperature()
        else:
            state = simulation.context.getState(getEnergy=True)
            temperature = (
                2 * state.getKineticEnergy() / (dof * MOLAR_GAS_CONSTANT_R)
                )
        tempreature_k = temperature.value_in_unit(kelvin)
        return tempreature_k

    @staticmethod
    def remove_water_and_ions(input_filepath: str) -> str:
        """Remove water from the output of an explicit solvent run.

        Uses the pdb-tools.pdb_delhetatm() module to do this

        NOTE: May produce an issue with TER atoms (Water + ions)
        
        NOTE2: May lead to issue with modified AA e.g: phsophoserine, etc...

        Parameters
        ----------
        input_filepath : str
            Path to the file to be modified

        Return
        ----------
        output_filepath : str
            Path to the modified file not containing anymore HETATM
        """
        log.info(f"Removing water and ions from {input_filepath}")
        output_filepath = input_filepath.replace('.pdb', '_nosolv.pdb')
        with open(input_filepath, 'r') as fh:
            pdbYieldedLines = pdb_delhetatm.run(fh)
            with open(output_filepath, 'w') as writefh:
                writefh.writelines(pdbYieldedLines)
        return output_filepath

    def run_openmm(
            self,
            inputPDBfile: str,
            output_directory: str,
            solvent_model: str,
            replica_index: int,
            ) -> dict:
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

        Return
        ----------
        output_files : dict
            Dict holding paths to configurations generated along the simulation
        """
        # Initiate output dictionary
        output_files: dict[str, Union[str, list[str]]] = {}
        # Split output path
        fn, ext = os.path.splitext(self.output_filename)

        # Make sure something has to be done...
        if self.params["simulation_timesteps"] == 0:
            output_filepath = os.path.join(
                output_directory,
                f"{fn}_{replica_index}{ext}"
                )
            os.rename(inputPDBfile, output_filepath)
            output_files["final"] = output_filepath
            return output_files

        # Load pdb
        pdb = openmmpdbfile(inputPDBfile)
        forcefield = ForceField(self.params['forcefield'], solvent_model)
        # system setup
        nonbondedMethod = NoCutoff if self.params["implicit_solvent"] else PME
        # should give an ERROR when system_constraints = 'None'.
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=nonbondedMethod,
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

        # Add the restrain force to the system
        # self._gen_centroid_forces(pdb)

        # setting integrator seed
        replica_seed = self.params["iniseed"] + replica_index
        if replica_seed != 0:
            log.info(f"setting non-random seed = {replica_seed}")
            integrator.setRandomNumberSeed(replica_seed)
        # retrieving integrator seed
        get_seed = integrator.getRandomNumberSeed()
        log.info(f"simulation seed {get_seed} for {self.model.file_name}")

        # simulation
        simulation = Simulation(pdb.topology, system, integrator)
        statereporter = statedatareporter.StateDataReporter(
            f'{self.directory_dict["simulation_stats"]}/observables_{self.identificator}_{replica_index}.dat', # noqa : E501
            50,
            step=True,
            totalEnergy=True,
            potentialEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            )
        simulation.reporters.append(statereporter)

        platform = simulation.context.getPlatform()
        log.info(f"Running on platform {platform.getName()}")
        log.info(f"Estimated platform speed: {platform.getSpeed()}")
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        log.info(f"Running NVT equilibration for {self.model.file_name}")

        # warming up the system with a NVT equilibration (canonical ensemble)
        eq_steps = self.params["equilibration_timesteps"]
        log.info(
            f"Warming up the system to {self.params['temperature_kelvin']}"
            f" K in ~{eq_steps} steps"
            )
        # In how many iterations are we inceasing the temperature ?
        nvt_sims = 10
        # How much should be increase the temperature at each iteration
        delta_temp = self.params["temperature_kelvin"] / nvt_sims
        integrator.setFriction(100 / picosecond)
        for n in range(1, nvt_sims + 1):
            qtemp = n * delta_temp
            integrator.setTemperature(qtemp * kelvin)
            simulation.step(int(eq_steps / nvt_sims))
        # Set final temperature and temperature coupling parameters
        integrator.setTemperature(self.params["temperature_kelvin"] * kelvin)
        integrator.setFriction(1 / picosecond)
        # Try to obtain the degree of freedom
        try:
            dof = statereporter._dof
        except AttributeError:
            # may not exists if in implicit solvent..
            dof = system.getNumParticles() * 3
        # Makes sure temperature of the system is reached
        self._stabilize_temperature(
            simulation,
            self.params["temperature_kelvin"],
            dof,
            tolerance=5.0,
            steps=50
            )
        log.info(
            f"Temperature of {self.params['temperature_kelvin']} K "
            f"reached in {simulation.context.getStepCount()} steps"
            )

        # Write
        eq_pdb_filepath = os.path.join(
            self.directory_dict["intermediates"],
            f"{fn}_eq_{replica_index}{ext}"
            )
        self._write_current_pdb(simulation, eq_pdb_filepath)
        output_files["equilibrated"] = eq_pdb_filepath
        
        # NPT simulation (isothermal-isobaric ensemble)
        _barostat = MonteCarloBarostat(  # noqa : F841
            1 * atmosphere,
            self.params["temperature_kelvin"] * kelvin
            )
        # FIXME : Once the MonteCarloBarostat will no more split chains
        # appart in different periodic boxes, please uncomment next line
        # simulation.system.addForce(_barostat)

        # Reinitialize the simulation
        simulation.context.reinitialize(True)

        # Running real simulation
        log.info(
            f"Running simulation for {self.model.file_name} for "
            f"{self.params['simulation_timesteps']} timesteps"
            )

        # Save intermediates
        if self.params["save_intermediate"]:
            # Initialize variables
            done_steps = 0
            output_files["intermediates"] = []
            steps_intermediate = int(
                self.params["simulation_timesteps"]
                // self.params["save_intermediate"]
                )
            starting_step_inter = int(steps_intermediate / 2)

            # Do few steps
            simulation.step(starting_step_inter)
            done_steps += starting_step_inter
            # Define intermediates output file
            int_filepath = os.path.join(
                self.directory_dict["intermediates"],
                f"{fn}_{replica_index}_{done_steps}steps{ext}"
                )
            self._write_current_pdb(simulation, int_filepath)
            output_files["intermediates"].append(int_filepath)

            # Loop over intermediate steps
            for _i in range(self.params["save_intermediate"] - 1):
                simulation.step(steps_intermediate)
                done_steps += steps_intermediate
                # Define intermediates output file
                int_filepath = os.path.join(
                    self.directory_dict["intermediates"],
                    f"{fn}_{replica_index}_{done_steps}steps{ext}"
                    )
                self._write_current_pdb(simulation, int_filepath)
                output_files["intermediates"].append(int_filepath)

            # Do the remaining simulation steps
            simulation.step(self.params["simulation_timesteps"] - done_steps)
        else:
            simulation.step(self.params["simulation_timesteps"])

        # Write final pdb file
        output_filepath = os.path.join(
            output_directory,
            f"{fn}_{replica_index}{ext}"
            )
        self._write_current_pdb(simulation, output_filepath)
        output_files["final"] = output_filepath
        # Log end of OpenMM simulation info
        log.info(f"OpenMM simulation successful for: {output_filepath}")
        # Return data
        return output_files

    @staticmethod
    def _write_current_pdb(
            simulation: Simulation,
            output_path: str,
            reference_pdb: Optional[list] = None,
            ) -> None:
        """Write the current position of atoms in the simulation.

        Parameters
        ----------
        simulation: :py:class:`openmm.app.simulation.Simulation`
            An openmm simulation class

        output_path : str
            Path of the file where to write the pdb file
        """
        # Retrieve last state data
        current_state = simulation.context.getState(
            getPositions=True,
            # enforcePeriodicBox=True
            )
        # Write pdb file
        openmmpdbfile.writeFile(
            simulation.topology,
            current_state.getPositions(),
            open(output_path, "w"),
            keepIds=True,
            )

    def generate_output_ensemble(
            self,
            replicas_outputs: list[dict[str, Union[str, list[str]]]],
            ) -> str:
        """Combine all replicas outputs into a single file.

        Parameters
        ----------
        replicas_outputs: list[dict[str, Union[str, list[str]]]]
            List of dictionnary mapping to all simulation replicas outputs.

        Returns
        -------
        final_openmm_output: str
            Path to the ensemble (composed of the equilibrated, intermediates,
            and final structure), for all replicas, generated by the module.
        """
        # Combines all outputs in a list
        all_files: list[str] = []
        for sample in replicas_outputs:
            for _outputtype, outputpath in sample.items():
                if isinstance(outputpath, str):
                    all_files.append(outputpath)
                elif isinstance(outputpath, list):
                    all_files += outputpath
        # Generate ensemble
        ensemble_filepath = os.path.join(
            self.directory_dict["openmm_output"],
            self.output_filename,
            )
        log.info(
            f"Generating ensemble from {len(all_files)} samples: "
            f"{ensemble_filepath}"
            )
        ensemble = pdb_mkensemble.run(all_files)
        with open(ensemble_filepath, "w") as wfile:
            for line in ensemble:
                wfile.write(line)
        # Remove single structure files
        for non_ensemble_fp in all_files:
            os.remove(non_ensemble_fp)
        return ensemble_filepath

    def run(self) -> None:
        """Run openmm simulation."""
        # using pdbfixer
        self.openmm_pdbfixer()
        # explicit solvent simulation
        if not self.params["implicit_solvent"]:
            # define solvent model and create solvation box
            solvent_model = self.params["explicit_solvent_model"]
            # create solvation box
            try:
                pdbPath = self.create_solvation_box(solvent_model)
            except ModuleError as _err:
                return
            # Eauilibrate this box around system
            pdbPath = self.equilibrate_solvation_box(pdbPath, solvent_model)
            # set variables
            output_folder = self.directory_dict["md_raw_output"]
            solvent = "explicit solvent"
        # implicit solvent simulation
        else:
            pdbPath = os.path.join(
                self.directory_dict["pdbfixer"],
                self.model.file_name
                )
            output_folder = self.directory_dict["openmm_output"]
            solvent_model = self.params["implicit_solvent_model"]
            solvent = "implicit solvent"

        log.info(f"starting {solvent} openMM simulation with file: {pdbPath}")
        
        # Loop over samplings/replicas
        replicas_outputs: list[dict[str, Union[str, list[str]]]] = []
        for replica_ind in range(self.params["sampling_factor"]):
            replica_struct = self.run_openmm(
                pdbPath,
                output_folder,
                solvent_model,
                replica_ind + 1,
                )
            replicas_outputs.append(replica_struct)

        # Solvent removal procedure
        if not self.params["implicit_solvent"] and \
                not self.params["keep_solvent"]:
            # Loop over replicas
            for replica_structs in replicas_outputs:
                for outputtype in replica_structs.keys():
                    # Remove solvent molecules from outputs
                    if isinstance(replica_structs[outputtype], list):
                        nosolv_out = [
                            self.remove_water_and_ions(pdb)
                            for pdb in replica_structs[outputtype]
                            ]  # type: ignore
                    else:
                        nosolv_out = self.remove_water_and_ions(
                            replica_structs[outputtype]
                            )  # type: ignore
                    # Modify paths in files mapper
                    replica_structs[outputtype] = nosolv_out

        # Post-process results and place them in appropriate directory
        if self.params["generate_ensemble"]:
            # Generate ensemble
            self.generate_output_ensemble(replicas_outputs)
        else:
            # Generate list of final structures only
            for replica_mapper in replicas_outputs:
                raw_output = replica_mapper["final"]
                final_output = raw_output.replace(
                    self.directory_dict["md_raw_output"],
                    self.directory_dict["openmm_output"]
                    )
                os.rename(raw_output, final_output)
