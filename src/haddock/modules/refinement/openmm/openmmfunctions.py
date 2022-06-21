from openmm import *
from openmm.app.pdbfile import PDBFile as openmmpdbfile
from openmm.app import *
from openmm.unit import *
from pdbfixer import PDBFixer
from pdbtools import pdb_delresname

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

def Does_pdb_contain_xray_crystallography_cell_data(pdb_path):
    openmmpdb = openmmpdbfile(pdb_path)
    modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
    if(modeller.topology.getUnitCellDimensions() is None):
        return False
    else:
        return True

# Order of PDBFixer method calls is important. For further info see: https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html
# Solvent is added using the OpenMM modeller in the createSolvationBox function.
def OpenmmPdbfixer(haddockmodule, pdb_path, pdb_filename, pdbfixer_result_directory):
    haddockmodule.log(f'Fixing pdb: {pdb_filename}')
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # fixer.addMissingHydrogens(7.0)
    openmmpdbfile.writeFile(fixer.topology, fixer.positions, open(os.path.join(pdbfixer_result_directory, pdb_filename), 'w'))

def createSolvationBox(haddockmodule, pdb_FilePath, pdb_file_name, buildPdbPath, forcefield, waterModel, contains_xray_crystallography_cell_data):
    try:
        if('.cif' in pdb_FilePath):
            openmmpdb = PDBxFile(pdb_FilePath)
        else:
            openmmpdb = openmmpdbfile(pdb_FilePath)
    
        modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
        usedForcefield = ForceField(forcefield, waterModel)
        modeller.deleteWater()
        modeller.addHydrogens(usedForcefield, pH=7.0)
      
        if contains_xray_crystallography_cell_data:
            modeller.addSolvent(usedForcefield, padding=haddockmodule.params['contains_xray_crystallography_cell_data_solvationbox_padding_nm']*nanometers)
        else:
            boxsizeIn_nm = haddockmodule.params['solvent_boxsize_nm']
            modeller.addSolvent(usedForcefield, boxSize=Vec3(boxsizeIn_nm, boxsizeIn_nm, boxsizeIn_nm)*nanometers, neutralize=False)

        # Add required extra particles for forcefield, e.g. Drude particles.
        if(haddockmodule.params['add_extra_particles_for_forcefield']):
            modeller.addExtraParticles(forcefield)

        if(haddockmodule.params['move_atoms_to_solvationbox_center']):
            AtomPositions = moveAtoms(modeller.positions, 10)
        else:
            AtomPositions = modeller.positions
        openmmpdbfile.writeFile(modeller.topology, AtomPositions, open(os.path.join(buildPdbPath, pdb_file_name), 'w'))

    except Exception as e:
        haddockmodule.log(f"BUILDING SOLVATION BOX ERROR: {e}")

def removeWaterMolecules(filePath, filename, destinationDirectory):
    with open(filePath, 'r') as fh:
        pdbYieldedLines = pdb_delresname.run(fh, 'HOH')
        with open(os.path.join(destinationDirectory, filename), 'w') as writefh:
            writefh.writelines(pdbYieldedLines)


def runOpenMM(haddockmodule, inputPDBfile, pdbFilename, outputDirectory, intermediate_structure_directory, forcefield_model, solvent_model):
    pdb = openmmpdbfile(inputPDBfile)
    forcefield = ForceField(forcefield_model, solvent_model)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=globals()[haddockmodule.params['constraints']], removeCMMotion=haddockmodule.params['remove_center_of_mass_motion'], rigidWater=haddockmodule.params['rigid_water']) #ERROR when system_constraints = 'None'.
    integrator = LangevinMiddleIntegrator(haddockmodule.params['temperature_kelvin']*kelvin, 1/picosecond, haddockmodule.params['timestep_ps']*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    # haddockmodule.log(f'Seed used by simulation integrator: {simulation.integrator.getRandomNumberSeed()}')
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    simulation.step(haddockmodule.params['equilibration_timesteps'])
    simulation.reporters.append(PDBReporter(os.path.join(outputDirectory, pdbFilename), haddockmodule.params['simulation_timesteps']))
    if(haddockmodule.params['save_intermediate_simulation_structures']):
        simulation.reporters.append(PDBReporter(os.path.join(intermediate_structure_directory, pdbFilename), haddockmodule.params['save_intermediate_simulation_structures_after_number_of_timesteps']))
    # simulation.reporters.append(StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True, temperature=True)) # Report system state of the simulation.
    simulation.step(haddockmodule.params['simulation_timesteps'])
    haddockmodule.log(f'OpenMM simulation succesful for: {os.path.join(outputDirectory, pdbFilename)}')