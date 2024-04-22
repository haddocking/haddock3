"""haddock3-restraints passive_from_active subcommand.

Given a list of active_residues and a PDB structure, it will return a list of
surface exposed passive residues within a 6.5A radius from the active residues.

When provided with a list of surface residues, it will filter the list for those
that are within 6.5A from the active residues.

Usage:
    haddock3-restraints passive_from_active <pdb_file> <active_list> 
        [-c <chain_id>]
        [-s <surface_list>]
        [--cutoff <relative_accessibility_cutoff>]
        [-n <distance_defining_neighboring>]
"""
from pathlib import Path

import numpy as np
import sys
from Bio.PDB import PDBParser, NeighborSearch

from haddock.core.typing import Optional
from haddock.libs.librestraints import get_surface_resids


def add_pass_from_act_arguments(pass_from_act_subcommand):
    """Add arguments to the pass_from_act subcommand."""
    pass_from_act_subcommand.add_argument(
        "structure",
        type=str,
        help="input PDB structure.",
        )

    pass_from_act_subcommand.add_argument(
        "active_list",
        help="List of active residues IDs (int) separated by commas",
        type=str,
        )

    pass_from_act_subcommand.add_argument(
        "-c",
        "--chain-id",
        help="Chain id to be used in the PDB file (default: All)",
        required=False,
        )

    pass_from_act_subcommand.add_argument(
        "-s",
        "--surface-list",
        help="List of surface residues IDs (int) separated by commas",
        required=False,
        type=str,
        )

    pass_from_act_subcommand.add_argument(
        "--cutoff",
        help="Relative cutoff for sidechain accessibility",
        required=False,
        default=0.15,
        type=float,
        )
    
    pass_from_act_subcommand.add_argument(
        "-n",
        "--neighbors-dist",
        help="Distance defining a neighbor from an active residue",
        required=False,
        default=6.5,
        type=float,
        )

    return pass_from_act_subcommand


def passive_from_active(
        structure,
        active_list,
        chain_id: Optional[str] = None,
        surface_list: Optional[str] = None,
        cutoff: float = 0.15,
        neighbors_dist: float = 6.5,
        ):
    """Get the passive residues."""

    # Parse the PDB file
    if Path(structure).exists():
        try:
            p = PDBParser(QUIET=True)
            s = p.get_structure('pdb', structure)
        except Exception as e:
            print(f'Error while parsing the PDB file: {e}')
            sys.exit(1)
    else:
        print(f'File not found: {structure}')
        sys.exit(1)

    try:
        if chain_id:
            atom_list = [a for a in s[0][chain_id].get_atoms()]
        else:
            atom_list = [a for a in s[0].get_atoms()]
    except KeyError as e:
        print(
            f'Chain {chain_id} does not exist in the PDB file {structure}'
            ', please enter a proper chain id'
            )
        sys.exit(1)

    try:
        active_list = [int(res) for res in active_list.split(',')]
        act_atoms = [
            a.get_coord()
            for a in atom_list
            if a.parent.id[1] in active_list
            ]
    except:
        print(
            'The list of active residues must '
            'be provided as a comma-separated list of integers'
            )
        sys.exit(1)

    try:
        if surface_list:
            surface_resids = [int(res) for res in surface_list.split(',')]
        else:
            surface_resids = get_surface_resids(s, cutoff=cutoff * 100)
    except Exception as e:
        print(f"There was an error while calculating surface residues: {e}")
        sys.exit(1)

    # Search for Neighbors
    ns = NeighborSearch(atom_list)
    neighbors = []
    for a in act_atoms:
        neighbors.append(
            ns.search(a, neighbors_dist, "R")  # HADDOCK used 6.5A as default
            )

    passive_list = set()
    for n in neighbors:
        for r in n:
            passive_list.add(r.id[1])
    tmp = passive_list & set(surface_resids)
    passive_list = tmp - set(active_list)
    print(' '.join([str(r) for r in sorted(passive_list)]))
    return
