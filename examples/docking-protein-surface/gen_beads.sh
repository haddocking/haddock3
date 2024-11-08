#!/bin/bash
#========================================================================================================================#
# Helper script to generate z-restraints for the example docking-protein-surface.
#
#  The necessary files are already provided in the `data/` directory. The script is provided as an example of how to
#   generate the z-restraints, in case you want to generate them for other proteins.
#
# Usage: ./gen_beads.sh
#
#========================================================================================================================#

# Generates only one plan based on one selection of residues
haddock3-restraints z_surface_restraints --pdb data/protein.pdb --residues 19,83,145,119,167 --output data/one_z_surface_selection

# Generates two plans based on the two residues selections
haddock3-restraints z_surface_restraints --pdb data/protein.pdb --residues 19,83,145,119,167 98,101,126,129 --output data/protein_z_restraints

# Generates two plans based on the two residues selections and reduce both area of the surface and space between beads
haddock3-restraints z_surface_restraints --pdb data/protein.pdb --residues 19,83,145,119,167 98,101,126,129 --output data/protein_z_restraints_smallspacing --spacing 4 --x-size 20 --y-size 20

# Generates three plans based on the three residues selections
haddock3-restraints z_surface_restraints --pdb data/protein.pdb --residues 19,83,145,119,167 98,101,126,129 23,62,87,111,116,153,163 --output data/three_z_restraints

#========================================================================================================================#