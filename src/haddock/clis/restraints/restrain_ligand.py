"""haddock3-restraints restrain_ligand subcommand.

Given an input PDB file and a residue name (ligands), the tool will create 
unambiguous restraints to keep this ligand in place during the refinements.


Usage:
    haddock3-restraints restrain_ligand <pdb_file> <ligand_name> -r <radius> -d <deviation> -n <max_nb_restraints>

positional arguments:
  pdb_file              Input PDB file.
  ligand_name           Residue name.

options:
  -h, --help            show this help message and exit
  -r RADIUS, --radius RADIUS
                        Radius used for neighbors search around center of mass of ligand.
                        (default: 10.0)
  -n MAX_RESTRAINTS, --max-restraints MAX_RESTRAINTS
                        Maximum number of restraints to return. (default: 200)
  -d DEVIATION, --deviation DEVIATION
                        Allowed deviation from actual distance. (default: 1.0)
"""

import os
import sys
import random
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch

from haddock.core.typing import Union


def add_restrain_ligand_arguments(restrain_ligand_subcommand):
    """Add arguments to the random removal subcommand."""
    restrain_ligand_subcommand.add_argument(
        "pdb_file",
        type=str,
        help="Input PDB file.",
        )
    restrain_ligand_subcommand.add_argument(
        "ligand_name",
        type=str,
        help="Name of the residue for which restraints must be generated.",
        )
    restrain_ligand_subcommand.add_argument(
        "-r",
        "--radius",
        help=(
            "Radius, in Angstrom, used for neighbors search around "
            "center of mass of ligand. (default: %(default)s)"
            ),
        required=False,
        default=10.0,
        type=float,
        )
    restrain_ligand_subcommand.add_argument(
        "-d",
        "--deviation",
        help=(
            "Allowed deviation from actual distance, in Angstrom. "
            "(default: %(default)s)"
            ),
        required=False,
        default=1.0,
        type=float,
        )
    restrain_ligand_subcommand.add_argument(
        "-n",
        "--max-restraints",
        help="Maximum number of restraints to return. (default: %(default)s)",
        required=False,
        default=200,
        type=int,
        )
    return restrain_ligand_subcommand


def restrain_ligand(
        pdbfile: Union[str, Path],
        ligand_name: str,
        radius: float = 10.0,
        deviation: float = 1.0,
        max_restraints: int = 20,
        ) -> str:
    """Generate unambiguous restraints to keep a ligand in place.

    Parameters
    ----------
    pdbfile : Union[str, Path]
        Path to a PDB file containing a ligand
    ligand_name : str
        Residue name of the ligand to work with
    radius : float, optional
        Radius used for neighbors search around center of mass of ligand, by default 10.0
    deviation : float, optional
        Allowed deviation from actual distance, by default 1.0
    max_restraints : int, optional
        Maximum number of restraints to return, by default 20

    Returns
    -------
    unambig_str : str
        The actual unambiguous restraints as a string.
    """
    # Read in structure
    pdb_parser = PDBParser(QUIET=1)
    structure = pdb_parser.get_structure("", pdbfile)

    # Remove hydrogens
    for atom in structure.get_atoms():
        if atom.element == "H":
            res = atom.parent
            res.detach_child(atom.name)

    # Try to find the ligand in the structure
    ligand = None
    for residue in structure.get_residues():
        if residue.resname.strip() == ligand_name.strip():
            ligand = residue
            break
    # Case where the ligand is not found in the structure
    if not ligand:
        print(f"[!!] Ligand residue '{ligand_name}' not found in structure")
        sys.exit(1)

    # Calculate center of mass of the ligand
    ligand_com = list(map(lambda x: sum(x)/len(x), zip(*[at.coord for at in ligand])))
    ligand_com = np.asarray(ligand_com, dtype=np.float32)

    # Create a selection of aminoacid/nucleotide atoms
    # (excl. waters, other ligands, etc)
    # also filters atoms that are from the queried ligand
    sel_atoms = [
        at for at in structure.get_atoms()
        if at.parent.id[0] == ' ' and at.parent != ligand
        ]
    # Perfom neighbor search on this selection
    ns = NeighborSearch(sel_atoms)
    neighbors = ns.search(ligand_com, radius, level="R")  # 10A radius, return residues

    # Calculate residue closer to each ligand atom and the respective distance
    ligand_atoms = ligand.child_list
    min_dist_list, _seen = [], set()

    # Loop over ligand atoms
    for l_atm in ligand_atoms:
        distances = []
        # Loop over neighbors residues and atoms
        for residue_atoms in neighbors:
            for r_atm in residue_atoms:
                # Compute distance and hold it
                distances.append((r_atm, l_atm, r_atm - l_atm))
        # Sort list by distances
        distances.sort(key=lambda x: x[-1])
        # Loop over sorted distances
        for closest_candidate in distances:
            candidate_residue = closest_candidate[0].parent
            # One restraint per residue to keep the number of restraints small
            # If a residue is already used, take next one
            if candidate_residue not in _seen:
                min_dist_list.append(closest_candidate)
                _seen.add(candidate_residue)
                break

    # Output
    assign_str_template = (
        "assign (segi {receptor_chainid:4s} and resi {receptor_resid:4d} "
        "and name {receptor_atname:6s}){linesep:s}"
        "     (segi {ligand_chainid:4s} and resi {ligand_resid:4d} "
        "and name {ligand_atname:6s}) "
        "{distance:6.3f} {deviation:.2f} {deviation:.2f}{linesep:s}"
        )

    _unambig_str_list: list[str] = []
    # Loop over all min distances
    for dist in min_dist_list:
        r_at, l_at, d = dist
        # Build assign statement
        assign_str = assign_str_template.format(
            receptor_chainid=r_at.parent.parent.id,
            receptor_resid=r_at.parent.id[1],
            receptor_atname='"' + r_at.name + '"',
            ligand_chainid=l_at.parent.parent.id,
            ligand_resid=l_at.parent.id[1],
            ligand_atname='"' + l_at.name + '"',
            distance=d,
            deviation=deviation,
            linesep=os.linesep,
            )
        _unambig_str_list.append(assign_str)
    
    # Limit the number of restraints
    if max_restraints < len(_unambig_str_list):
        random.seed(420)
        unambig_str_list = random.sample(_unambig_str_list, max_restraints)
    else:
        unambig_str_list = _unambig_str_list
    
    # Concatenate into a single string
    unambig_str = f"! Restraints to fix {ligand_name} in its initial position{os.linesep}"
    unambig_str += "".join(unambig_str_list)
    return unambig_str


def main(
        pdb_file: Union[str, Path],
        ligand_name: str,
        radius: float = 10.0,
        deviation: float = 1.0,
        max_restraints: int = 20,
        ) -> None:
    """Simple wrapper of the restrain_ligand function."""
    unambig_str = restrain_ligand(
        pdb_file, ligand_name,
        radius=radius,
        deviation=deviation,
        max_restraints=max_restraints,
        )
    print(unambig_str)
