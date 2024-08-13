# Generates only one plan based on one selection of residues
./haddock-restraints z \
  --residues 19,83,145,119,167 \
  data/protein.pdb \
  data/one_z_surface_selection_new.pdb \
  20 \
  6 > data/one_z_surface_selection_beads_new.tbl

# Generates two plans based on the two residues selections
./haddock-restraints z \
  --residues 19,83,145,119,167 \
  --residues 98,101,126,129 \
  data/protein.pdb \
  data/protein_z_restraints_beads_new.pdb \
  20 \
  6  > data/protein_z_restraints_new.tbl

# Generates two plans based on the two residues selections and reduce both area of the surface and space between beads
./haddock-restraints z \
  --residues 19,83,145,119,167 \
  --residues 98,101,126,129 \
  data/protein.pdb \
  data/protein_z_restraints_smallspacing_beads_new.pdb \
  4 \
  6  > data/protein_z_restraints_smallspacing_new.tbl

# Generates three plans based on the three residues selections
./haddock-restraints z \
  --residues 19,83,145,119,167 \
  --residues 98,101,126,129 \
  --residues 23,62,87,111,116,153,163 \
  data/protein.pdb \
  data/three_z_restraints_beads_new.pdb \
  20 \
  6  > data/three_z_restraints_new.tbl
