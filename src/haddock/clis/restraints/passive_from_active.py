"""haddock3-restraints passive_from_active subcommand.

Given a list of active_residues and a PDB structure, it will return a list of
surface exposed passive residues within a 6.5A radius from the active residues.

When provided with a list of surface residues, it will filter the list for those
that are within 6.5A from the active residues.

Usage:
    haddock3-restraints passive_from_active <pdb_file> <active_list> [-c <chain_id>] [-s <surface_list>]
"""

import sys

from haddock.libs.librestraints import passive_from_active_raw


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

    active = [int(res) for res in active_list.split(',')]
    surface = []
    if surface_list:
        surface = [int(res) for res in surface_list.split(',')]

    try:
        passive = passive_from_active_raw(
            structure, active, chain_id, surface
        )
    except Exception as e:
        print(e)
        sys.exit(1)

    print(' '.join([str(r) for r in passive]))
    return
