# Generates to plans based on the two residues selections
haddock3-restraints z_surface_restraints --pdb data/protein.pdb --residues 19,83,145,119,83,145,167 98,101,126,129 --output data/protein_z_restraints

# Generates to plans based on the two residues selections and reduce the space between beads
haddock3-restraints z_surface_restraints --pdb data/protein.pdb --residues 19,83,145,119,83,145,167 98,101,126,129 --output data/protein_z_restraints_smallspacing --spacing 4 --x-size 20 --y-size 20
