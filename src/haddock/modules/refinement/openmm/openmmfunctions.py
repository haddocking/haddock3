"""Operate function from OpenMM."""
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


def move_atoms(atomPositionList, nanometerBoxSize, addition=True):
    """Move atoms."""
    shiftedAtomPostions = []
    if addition:
        sign = 1.0
    else:
        sign = -1.0

    shift = sign * nanometerBoxSize / 2 * 10
    for Vector3 in atomPositionList:
        newAtomVec3PositionManual = Vec3(
            Vector3.x * 10 + shift,
            Vector3.y * 10 + shift,
            Vector3.z * 10 + shift,
            )
        shiftedAtomPostions.append(newAtomVec3PositionManual)

    return shiftedAtomPostions


def does_pdb_contain_xray_crystallography_cell_data(pdb_path):
    """Check of PDB contains x-ray cell data."""
    openmmpdb = openmmpdbfile(str(pdb_path))
    modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
    if modeller.topology.getUnitCellDimensions() is None:
        return False
    else:
        return True


# Order of PDBFixer method calls is important.
# For further info see: https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html  # noqa: E501
# Solvent is added using the OpenMM modeller in the createSolvationBox function.
def OpenmmPdbfixer(
        haddockmodule,
        pdb_path,
        pdb_filename,
        pdbfixer_result_directory,
        ):
    """Make OpenMM PDB Fixer."""
    haddockmodule.log(f'Fixing pdb: {pdb_filename}')
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # fixer.addMissingHydrogens(7.0)
    out_file = Path(pdbfixer_result_directory, pdb_filename)
    with open(out_file, 'w') as fout:
        openmmpdbfile.writeFile(fixer.topology, fixer.positions, fout)


def createSolvationBox(
        haddockmodule,
        pdb_FilePath,
        pdb_file_name,
        buildPdbPath,
        forcefield,
        waterModel,
        contains_xray_crystallography_cell_data
        ):
    """Create solvation box."""
    try:
        if pdb_FilePath.endswith('.cif'):
            openmmpdb = PDBxFile(pdb_FilePath)
        else:
            openmmpdb = openmmpdbfile(pdb_FilePath)

        modeller = Modeller(openmmpdb.topology, openmmpdb.positions)
        usedForcefield = ForceField(forcefield, waterModel)
        modeller.deleteWater()
        modeller.addHydrogens(usedForcefield, pH=7.0)

        if contains_xray_crystallography_cell_data:
            modeller.addSolvent(
                usedForcefield,
                padding=haddockmodule.params['contains_xray_crystallography_cell_data_solvationbox_padding_nm'] * nanometers,  # noqa: E501
                )
        else:
            boxsizeIn_nm = haddockmodule.params['solvent_boxsize_nm']
            modeller.addSolvent(
                usedForcefield,
                boxSize=Vec3(boxsizeIn_nm, boxsizeIn_nm, boxsizeIn_nm) * nanometers,  # noqa: E501
                neutralize=False,
                )

        # Add required extra particles for forcefield, e.g. Drude particles.
        if haddockmodule.params['add_extra_particles_for_forcefield']:
            modeller.addExtraParticles(forcefield)

        if haddockmodule.params['move_atoms_to_solvationbox_center']:
            AtomPositions = move_atoms(modeller.positions, 10)
        else:
            AtomPositions = modeller.positions

        output_file = Path(buildPdbPath, pdb_file_name)
        with open(output_file, 'w') as fout:
            openmmpdbfile.writeFile(
                modeller.topology,
                AtomPositions,
                fout,
                )

    except Exception as err:
        haddockmodule.log(f"BUILDING SOLVATION BOX ERROR: {err}")


def removeWaterMolecules(filePath, filename, destinationDirectory):
    """Remove water molecules."""
    with open(filePath, 'r') as fh:
        pdbYieldedLines = pdb_delresname.run(fh, 'HOH')
        with open(Path(destinationDirectory, filename), 'w') as writefh:
            writefh.writelines(pdbYieldedLines)


def runOpenMM(
        haddockmodule,
        inputPDBfile,
        pdbFilename,
        outputDirectory,
        intermediate_structure_directory,
        forcefield_model,
        solvent_model,
        ):
    """Run OpenMM."""
    pdb = openmmpdbfile(inputPDBfile)
    forcefield = ForceField(forcefield_model, solvent_model)

    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * nanometer,
        constraints=haddockmodule.params['constraints'],
        removeCMMotion=haddockmodule.params['remove_center_of_mass_motion'],
        rigidWater=haddockmodule.params['rigid_water'],
        )  # ERROR when system_constraints = 'None'.

    integrator = LangevinMiddleIntegrator(
        haddockmodule.params['temperature_kelvin'] * kelvin,
        1 / picosecond,
        haddockmodule.params['timestep_ps'] * picoseconds,
        )

    simulation = Simulation(pdb.topology, system, integrator)

    # haddockmodule.log(
    # 'Seed used by simulation integrator: '
    # f'{simulation.integrator.getRandomNumberSeed()}')
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    simulation.step(haddockmodule.params['equilibration_timesteps'])
    simulation.reporters.append(
        PDBReporter(
            str(Path(outputDirectory, pdbFilename)),
            haddockmodule.params['simulation_timesteps'],
            )
        )

    if haddockmodule.params['save_intermediate_simulation_structures']:
        simulation.reporters.append(
            PDBReporter(
                str(Path(intermediate_structure_directory, pdbFilename)),
                haddockmodule.params['save_intermediate_simulation_structures_after_number_of_timesteps'],  # noqa: E501
                )
            )

    # simulation.reporters.append(StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True, temperature=True)) # Report system state of the simulation.  # noqa: E501
    simulation.step(haddockmodule.params['simulation_timesteps'])
    haddockmodule.log(
        'OpenMM simulation succesful for: '
        f'{os.path.join(outputDirectory, pdbFilename)}'
        )
