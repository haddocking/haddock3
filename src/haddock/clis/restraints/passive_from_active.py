"""haddock3-restraints passive_from_active subcommand.

Given a list of active_residues and a PDB structure, it will return a list of
surface exposed passive residues within a 6.5A radius from the active residues.

When provided with a list of surface residues, it will filter the list for those
that are within 6.5A from the active residues.

Usage:
    haddock3-restraints passive_from_active <pdb_file> <active_list> [-c <chain_id>] [-s <surface_list>]
"""
from pathlib import Path

import numpy as np
import sys
from Bio.PDB import PDBParser, NeighborSearch

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

    return pass_from_act_subcommand


def passive_from_active(structure, active_list, chain_id=None, surface_list=""):
    """Get the passive residues."""

    # Parse the PDB file
    if Path(structure).exists():
        try:
            p = PDBParser(QUIET=True)
            s = p.get_structure('pdb', structure)
        except Exception as e:
            print('Error while parsing the PDB file: {0}'.format(e))
            sys.exit(1)
    else:
        print('File not found: {0}'.format(structure))
        sys.exit(1)

    try:
        if chain_id:
            atom_list = [a for a in s[0][chain_id].get_atoms()]
        else:
            atom_list = [a for a in s[0].get_atoms()]
    except KeyError as e:
        print('Chain {0} does not exist in the PDB file {1}, please enter a proper chain id'.
              format(chain_id, structure))
        sys.exit(1)

    try:
        active_list = [int(res) for res in active_list.split(',')]
        act_atoms = [a.get_coord() for a in atom_list if a.parent.id[1] in active_list]
    except:
        print('The list of active residues must be provided as a comma-separated list of integers')
        sys.exit(1)

    try:
        if surface_list:
            surface_list = [int(res) for res in surface_list.split(',')]
        else:
            surface_list = get_surface_resids(s)
    except Exception as e:
        print("There was an error while calculating surface residues: {}".format(e))
        sys.exit(1)

    ns = NeighborSearch(atom_list)
    neighbors = []
    for a in act_atoms:
        neighbors.append(ns.search(a, 6.5, "R"))  # HADDOCK used 6.5A as default

    passive_list = set()
    for n in neighbors:
        for r in n:
            passive_list.add(r.id[1])
    tmp = passive_list & set(surface_list)
    passive_list = tmp - set(active_list)
    print(' '.join([str(r) for r in sorted(passive_list)]))
    return 