"""haddock3-restraints restrain_bodies subcommand.

The restrain_bodies subcommand creates distance restraints to lock several 
chains together. Useful to avoid unnatural flexibility or movement due to
sequence/numbering gaps.

Usage::
	
	haddock3-restraints restrain_bodies <structure> [--exclude] [--verbose]
"""

import logging
import sys
import random

from haddock.libs.librestraints import read_structure, get_bodies, build_restraints, generate_tbl

# Set random seed to have reproducibility
random.seed(917)

def add_restrain_bodies_arguments(restraint_bodies_subcommand):
	restraint_bodies_subcommand.add_argument(
		"structure",
		help="The PDB structure to be restrained.",
		)
	
	restraint_bodies_subcommand.add_argument(
		"-e",
		"--exclude",
		help="Chains to exclude from the calculation.",
		required=False,
		type=str,
		)
	
	restraint_bodies_subcommand.add_argument(
		"-v",
		"--verbose",
		help="Tune verbosity of the output.",
		required=False,
		default=0,
		type=int,
		)
	
	restraint_bodies_subcommand


def restrain_bodies(structure, exclude=None, verbose=0):  # noqa: E501
	"""Create distance restraints to lock several chains together.
	
	Parameters
	----------
	structure : str
		The PDB structure to be restrained.

	exclude : str
		Chains to exclude from the calculation.
	
	verbose : int
		Tune verbosity of the output.
	"""
	if verbose== 1:
		logging.basicConfig(level=logging.INFO)
	elif verbose > 1:
		logging.basicConfig(level=logging.DEBUG)
	else:
		logging.basicConfig(level=logging.WARNING)

	# Main logic
	atom_lst = read_structure(structure, exclude=exclude)
	bodies = get_bodies(atom_lst)
	restraints = build_restraints(bodies)
	generate_tbl(atom_lst, restraints)
