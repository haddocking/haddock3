==================================================
    Peptide cyclisation protocol with HADDOCK3

This example workflow will take a peptide
and generate cyclised conformations in a two step
process, first using distance restraints to bring the
termini together, then rebuilding the topology to
create the covalent cyclic bond and refining again.
50 clusters are generated and the best model of
each is selected.

Protocol described in: [https://doi.org/10.1021/acs.jctc.2c00075](https://doi.org/10.1021/acs.jctc.2c00075)

==================================================

Two cyclic peptide examples are provided:

1) 1SFI, a 14 residue cyclic peptide with both backbone and disulphide bridge cyclisation

2) 3WNE, a 6 residue backbone cyclic peptide
