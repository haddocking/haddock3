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
./../../utils/haddock-restraints z \
  --residues 19,83,145,119,167 \
  data/protein.pdb \
  data/01_beads.pdb \
  20 \
  6 > data/01_z.tbl

# Generates two plans based on the two residues selections
./../../utils/haddock-restraints z \
  --residues 19,83,145,119,167 \
  --residues 98,101,126,129 \
  data/protein.pdb \
  data/02_beads.pdb \
  20 \
  6  > data/02_z.tbl

# Generates two plans based on the two residues selections and reduce both area of the surface and space between beads
./../../utils/haddock-restraints z \
  --residues 19,83,145,119,167 \
  --residues 98,101,126,129 \
  data/protein.pdb \
  data/02_beads_smallspacing.pdb \
  4 \
  6  > data/02_z_smallspacing.tbl

# Generates three plans based on the three residues selections
./../../utils/haddock-restraints z \
  --residues 19,83,145,119,167 \
  --residues 98,101,126,129 \
  --residues 23,62,87,111,116,153,163 \
  data/protein.pdb \
  data/03_beads.pdb \
  20 \
  6  > data/03_z.tbl

#========================================================================================================================#