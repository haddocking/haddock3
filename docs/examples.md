# Workflow Examples

The HADDOCK3 `examples/` directory contains various subdirectories and config files
corresponding to different types of complexes, scenarios and data.

1. [docking-antibody-antigen](#docking-antibody-antigen)
1. [docking-nanobody-antigen](#docking-nanobody-antigen)
1. [docking-protein-DNA](#docking-protein-dna)
1. [docking-protein-glycan](#docking-protein-glycan)
1. [docking-protein-homotrimer](#docking-protein-homotrimer)
1. [docking-protein-ligand](#docking-protein-ligand)
1. [docking-protein-ligand-shape](#docking-protein-ligand-shape)
1. [docking-protein-peptide](#docking-protein-peptide)
1. [docking-protein-protein](#docking-protein-protein)
1. [refine-complex](#refine-complex)
1. [scoring](#scoring)
1. [docking-multiple-ambig](#docking-multiple-ambig)
1. [general analysis tasks](#general-analysis-tasks)

Each directory contains both:

- __test__ runs config files (to test the various workflows, make sure the installation works and any changes to the code/scripts have not broken the machinery). Those are set up to run locally.
- __full__ runs config files with recommended parameter settings. Those runs are set up to be executed in "batch" mode using slurm (the `...full.cfg`) files. Examples making use of MPI are also provided in some cases, together with an associated job file that should be submitted to the slurm batch system (`...full-mpi.cfg` and `...full-mpi.job`). Make sure to adapt the full config files to your own system.

The following examples are currently provided:


## docking-antibody-antigen

An antibody-antigen docking example making use only of the knowledge of the hypervariables (HV) loops on the antibody to guide the docking. This is the same complex used in our [HADDOCK2.4 webserver tutorial](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-antibody-antigen/); refer to it for more details. Three different ways of using the knowledge of the HV loop residues are illustrated:

- __CDR-accessible__: using ambiguous distance restraints (AIRs) between the HV loop residues and the solvent-accessible residues of the antigen. Those are defined in the `ambig.tbl` file provided in the `data` directory
- __ranairCDR__: using random AIRs sampled from the HV loops on the antibody and the solvent-accessible residues on the antigen. The random AIRs are used during the rigidbody docking sampling phase and replaced by contact AIRs during refinement.
- __CDR-NMR-CSP__: using ambiguous distance restraints (AIRs) between the HV loop residues and the epitope residues identified by NMR.

Since antibodies consist of two separate chains, a few additional unambiguous restraints are defined between the antibody chains to keep them together during the flexible refinement stage. Those are defined in the `unambig.tbl` file in the `data` directory.

Three different protocols/workflows are illustrated:

1) 10000 rigidbody docking models, selection of top500 and flexible refinement + EM of those ([docking-antibody-antigen-CDR-accessible-full.cfg](../examples/docking-antibody-antigen/docking-antibody-antigen-CDR-accessible-full.cfg) and [`docking-antibody-antigen-ranairCDR-full.cfg`](../examples/docking-antibody-antigen/docking-antibody-antigen-ranairCDR-full.cfg))
2) 10000 rigidbody docking models, FCC clustering and selection of max 20 models per cluster followed by flexible refinement and EM ([docking-antibody-antigen-CDR-accessible-clt-full.cfg](../examples/docking-antibody-antigen/docking-antibody-antigen-CDR-accessible-clt-full.cfg) and [docking-antibody-antigen-ranairCDR-clt-full.cfg](../examples/docking-antibody-antigen/docking-antibody-antigen-ranairCDR-clt-full.cfg)`).
3) 1000 rigidbody docking models, selection of top200 and flexible refinement + EM of those (default sampling in the case in which the epitope has been determined by NMR: [docking-antibody-antigen-CDR-NMR-CSP-full.cfg](../examples/docking-antibody-antigen/docking-antibody-antigen-CDR-NMR-CSP-full.cfg))

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.

## docking-nanobody-antigen

A nanobody-antigen docking example making use of different levels of knowledge about the antigen epitope.

The standard HADDOCK workflow is used for all the cases, with 1000 `rigidbody` docking models, selection of top200 for flexible refinement (`flexref`) + energy minimisation (`emref`).

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.


## docking-protein-dna

A protein-DNA docking example making use of rather specific protein-DNA distance restraints defined in the `ambig.tbl` file in the `data` directory.

Three different protocols/workflows are illustrated:

1) 1000 rigidbody docking models, selection of top200 and flexible refinement + EM of those ([docking-protein-DNA-full.cfg](../examples/docking-protein-DNA/docking-protein-DNA-full.cfg))
2) 1000 rigidbody docking models, selection of top200 and flexible refinement + final refinement in explict solvent (water) of those ([docking-protein-DNA-mdref-full.cfg](../examples/docking-protein-DNA/docking-protein-DNA-mdref-full.cfg)
3) 1000 rigidbody docking models, FCC clustering and selection of max 20 models per cluster followed by flexible refinement and EM ([docking-protein-DNA-cltsel-full.cfg](../examples/examples/docking-protein-DNA/docking-protein-DNA-cltsel-full.cfg)).

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.

__Note__ the modified electrostatic treatment for DNA (`cdie` and `epsilon=78`) and the automatic definition of DNA restraints to maintain the double helix (`dnarest_on=true`).

## docking-protein-glycan

A protein-glycan docking example making use of the knowledge of the binding site on the protein to guide the docking. The conformation of the glycan has been obtained from the [GLYCAM webserver](http://glycam.org/), while the structure of the protein is taken from the PDB in its unbound form. In the proposed workflows, a clustering step is always performed after initial docking stage, so as to increase the diversity of the ensemble of models to be refined. Three different workflows are illustrated:

1) 1000 rigidbody docking models, RMSD clustering to select 50 clusters, flexible refinement of the top 5 models of each cluster, final RMSD clustering for cluster-based scoring ([docking-protein-glycan-full.cfg](../examples/docking-protein-glycan/docking-protein-glycan-full.cfg)). The RMSD clustering assumes a good knowledge of the interface, as the user has to define the residues involved in the binding site by means of the resdic_ parameter.
2) 1000 rigidbody docking models, interface-ligand-RMSD (`ilrmsd`) clustering to select 50 clusters, flexible refinement of the top 5 models of each cluster, final ilRMSD clustering for cluster-based scoring ([docking-protein-glycan-ilrmsd-full.cfg](../examples/docking-protein-glycan/docking-protein-glycan-ilrmsd-full.cfg)). The interface-ligand-RMSD clustering is a more general approach, as it does not require the user to define the residues involved in the binding site. The interface is automatically defined by the residues involved in the protein-glycan interaction in the input models.
3) 500 flexible docking runs + final RMSD clustering for cluster-based scoring [docking-flexref-protein-glycan-full.cfg](../examples/docking-protein-glycan/docking-flexref-protein-glycan-full.cfg). In this case, the rigidbody docking is skipped and the docking is performed at the flexible refinement level. In this case the flexible refinement has more steps than usual (`mdsteps_rigid = 5000`, `mdsteps_cool1 = 5000` and so on) and the glycan is defined as fully flexible (`fle_sta_1`, `fle_end_1`, `fle_seg_1`).

__Note__ the modified weight of the van der Waals energy term for the scoring of the rigidbody docking models (`w_vdw = 1.0`), as in the [protein-ligand example](#docking-protein-ligand).


## docking-protein-homotrimer

A protein homotrimer with C3 symmetry docking example making use of bioinformatics interface predictions from our [Whiscy webserver](https://wenmr.science.uu.nl/whiscy/) defined in the [1qu9_whiscy_air.tbl](../examples/docking-protein-homotrimer/data/1qu9_whiscy_air.tbl) file in the `data` directory. Note that this is a bound docking example since the protein monomers used for docking are taken from the crystal structure of the homotrimer. The workflow does however perform a flexible refinement of the models, which might lead to deviations from the bound conformation.

The [docking-protein-homotrimer-full.cfg](../examples/docking-protein-homotrimer/docking-protein-homotrimer-full.cfg) workflow consists of the generation of 2000 rigidbody docking models, selection of top200 and flexible refinement + EM of those. C3 and non-crystallographic (NCS) symmetry restraints are defined to guide the docking.

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.


## docking-protein-ligand

A protein-ligand docking example making use of the knowledge of the binding site on the protein to guide the docking.

As explained in our [protein-ligand HADDOCK2.4 tutorial](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-binding-sites/), in the rigidbody docking phase all residues of the binding site are defined as active to draw the ligand into it (the corresponding AIRs are defined in the [ambig-active-rigidbody.tbl](../examples/docking-protein-ligand/data/ambig-active-rigidbody.tbl) file in the `data` directory). For the flexible refinement only the ligand is defined as active and the binding site as passive to allow the ligand to explore the binding site (the corresponding AIRs are defined in the [ambig-passive.tbl](../examples/docking-protein-ligand/data/ambig-passive.tbl) file in the `data` directory).

The [docking-protein-ligand-full.cfg](../examples/docking-protein-ligand/docking-protein-ligand-full.cfg) workflow consists of the generation of 1000 rigidbody docking models, selection of top200 and flexible refinement of those.

__Note__ the modified weight of the van der Waals energy term for the scoring of the rigidbody docking models (`w_vdw = 1.0`) and the skipping of the high temperature first two stages of the simulated annealing protocol during the flexible refinement (`mdsteps_rigid = 0` and `mdsteps_cool1 = 0`). Parameter and topology files must be provided for the ligand (`ligand_param_fname = "data/ligand.param"` and `ligand_top_fname = "data/ligand.top"`). Those were obtained with a local version of PRODRG ([Schüttelkopf and van Aalten Acta Crystallogr. D 60, 1355−1363 (2004)](http://scripts.iucr.org/cgi-bin/paper?S0907444904011679)).

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.

## docking-protein-ligand-shape

A protein-ligand docking example making use of the knowledge of a template ligand (a ligand similar to the ligand we want to dock and bind to the same receptor). The template ligand information is used in the form of shape consisting of dummy beads and positioned within the binding site to which distance restraints are defined. More details about the method and the performance of the protocol when benchmarked on a fully unbound dataset
can be seen in our freely available [paper on JCIM](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c00796).

As explained in our [shape small molecule HADDOCK2.4 tutorial](https://www.bonvinlab.org/education/HADDOCK24/shape-small-molecule/), during the docking and refinement the protein and the shape are kept in their original positions (see the `mol_fix_origin_X` parameters in the config file) and ambiguous distance restraints between the ligand and the shape beads are defined. (The corresponding AIRs are defined in the `shape-restraints-from-shape-1.tbl` file in the `data` directory). This is effectively a three body docking. For the ligand an ensemble of 10 different conformations is provided as starting point for the docking (`ligand-ensemble.pdb` in the `data` directory). Please refer to our [shape small molecule tutorial](https://www.bonvinlab.org/education/HADDOCK24/shape-small-molecule/) for information on how to generate such an ensemble.

The [docking-protein-ligand-shape-full.cfg](../examples/docking-protein-ligand-shape/docking-protein-ligand-shape-full.cfg) workflow consists of the generation of 1000 rigidbody docking models with the protein and shape kept in their origin position, selection of top200 and flexible refinement of those.

__Note__ the modified weight of the van der Waals energy term for the scoring of the rigidbody docking models (`w_vdw = 1.0`). To allow the ligand to penetrate better into the binding site the intermolecular energy components are scaled down during the rigidbody docking phase (`inter_rigid = 0.001`). As for the protein-ligand example, parameter and topology files must be provided for the ligand (`ligand_param_fname = "data/ligand.param"` and `ligand_top_fname = "data/ligand.top"`). Those were obtained with a local version of PRODRG ([Schüttelkopf and van Aalten Acta Crystallogr. D 60, 1355−1363 (2004)](http://scripts.iucr.org/cgi-bin/paper?S0907444904011679)).


The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.


## docking-protein-peptide

The protein-peptide docking example makes use of the knowledge of the binding site on the protein to guide the docking. The active site residues are defined as active and the peptide as passive (the corresponding AIRs are defined in the `ambig.tbl` file in the `data` directory). This example follows the protocol described in our protein-peptide docking article ([Trellet et. al. PLoS ONE 8, e58769 (2013)](https://dx.plos.org/10.1371/journal.pone.0058769)). For the peptide, an ensemble of three conformations (alpha-helix, polyproline-II and extended) is provided as starting point for the docking. Those were built using PyMol (instructions on how to do that can be found [here](https://www.bonvinlab.org/education/molmod_online/simulation/#preparing-the-system)).

Three different workflows are illustrated:

1) 3000 rigidbody docking models, selection of top400 and flexible refinement + EM of those ([docking-protein-peptide-full.cfg](../examples/docking-protein-peptide/docking-protein-peptide-full.cfg)
2) 3000 rigidbody docking models, selection of top400 and flexible refinement + final refinement in explicit solvent (water) of those ([docking-protein-peptide-mdref-full.cfg](../examples/docking-protein-peptide/docking-protein-peptide-mdref-full.cfg)
3) 3000 rigidbody docking models, FCC clustering and selection of max 20 models per cluster followed by flexible refinement and EM ([docking-protein-peptide-cltsel-full.cfg](../examples/docking-protein-peptide/docking-protein-peptide-cltsel-full.cfg)).

__Note__ how the peptide is defined as fully flexible for the refinement phase (`fle_sta_1`, `fle_end_1`, `fle_seg_1`) and dihedral angle restraints are automatically defined to maintain secondary structure elements (`ssdihed = "alphabeta"`)

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.


## docking-protein-protein

The protein-protein docking example makes use of the NMR chemical shift perturbation data providing information on the residues of binding site to guide the docking. The NMR-identified residues are defined as active with their surface neighbors as passive (the corresponding AIRs are defined in the [e2a-hpr_air.tbl](../examples/docking-protein-protein/data/e2a-hpr_air.tbl) file in the `data` directory). This system is the same as described in our [HADDOCK2.4 basic protein-protein docking tutorial](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-protein-protein-basic/). For the second molecule (HPR), an ensemble of 10 conformations (taken from the NMR solution structure of this protein) is used as starting point for the docking. Refer to above tutorial for more details about the system and restraints.

Three different workflows are illustrated:

1) 1000 rigidbody docking models, selection of top200 and flexible refinement + EM of those ([docking-protein-protein-full.cfg](../examples/docking-protein-protein/docking-protein-protein-full.cfg))
2) 1000 rigidbody docking models, selection of top200 and flexible refinement + final refinement in explicit solvent (water) of those ([docking-protein-protein-mdref-full.cfg](../examples/docking-protein-protein/docking-protein-protein-mdref-full.cfg))
3) 1000 rigidbody docking models, FCC clustering and selection of max 20 models per cluster followed by flexible refinement and EM ([docking-protein-protein-cltsel-full.cfg](../examples/docking-protein-protein/docking-protein-protein-cltsel-full.cfg)).

In this example all parameters are left to their default settings, except for manually defining the histidines protonation states.

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.


## refine-complex

This example illustrates the refinement of a complex. In this case (workflow [refine-complex-test.cfg](../examples/refine-complex/refine-complex-test.cfg)) the molecules are kept in their original positions and the complex is subjected to a short flexible refinement in explicit solvent with the `mdref` module. The same complex as for the `docking-protein-protein` example is used. The molecules are defined separately in the config file (and could consist each of an ensemble, provided the two ensembles have exactly the same number of models).

In this example all parameters are left to their default settings, except for manually defining the histidines' protonation states and setting the `sampling_factor` to 10, which means that from each starting complex 10 models will be generated with different random seeds for initiating the molecular dynamics phase.

The `caprieval` module is called at the end to assess the quality of the models with respect to the known reference structure.


## scoring

This example illustrates the use of HADDOCK3 for scoring purposes. In contrast to HADDOCK2.4, HADDOCK3 can score a heterogenous set of complexes within one run/workflow. In this example, four different types of complexes are scored within the same workflow:

- an ensemble of 5 models taken from CAPRI Target161
- a protein-DNA complex (model taken from our protein-DNA docking example)
- two models of a protein-protein complex (taken from our protein-protein docking example)
- a homotrimer model (taken from our protein-homotrimer docking examples)

Two scoring workflows are illustrated:

1) Only a short energy minimisation is performed on each model ([emscoring-test.cfg](../examples/scoring/emscoring-test.cfg)).
2) A short molecular dynamics simulation in explicit solvent (water) is performed on each model. In that case contact AIRs (`contactairs = true`), dihedral angle restraints on secondary structure element (`ssdihed = alphabeta`) and DNA restraints (`dnarest_on = true`) are automatically defined ([mdscoring-test.cfg](../examples/scoring/mdscoring-test.cfg)).
3) An example scoring pipeline using in the CAPRI55 competition ([capri-scoring-test.cfg](../examples/scoring/capri-scoring-test.cfg)), where energy minimisation is followed by FCC clustering and selection of the top 2 models per cluster. Then a short molecular dynamics simulation in explicit solvent (water) is performed on each model and the models are clustered again.

The model listings with their associated HADDOCK scores can be found in a `.tsv` file in the stage 01 directory of the respective runs.

## docking-multiple-ambig

This example shows how to use HADDOCK3 when several restraint files are available. The example ([docking-multiple-tbls-test.cfg](../examples/docking-multiple-ambig/docking-multiple-tbls-test.cfg)) is built upon the results obtained running [arctic3d](https://github.com/haddocking/arctic3d) on two proteins forming the complex `2GAF`. The presence of multiple interfaces in both structures allows to define several `.tbl` ambiguous restraint files to be used in the calculations. At first, these files must be compressed in a `.tbl.tgz` archive. During the workflow, the HADDOCK3 machinery unzips the archive and evenly assigns each `.tbl` file to a number of models to be generated. Even if only one sixth of the restraint files contain reasonable information on the interface, HADDOCK3 is still able to retrieve good docking models in the best-scoring positions.

__Note__ how the information about restraint files is propagated during the workflow (`previous_ambig = true` for `flexref` and `emref` modules), so that each model is always refined with its corresponding `.tbl` file.

Importantly, in the [docking-multiple-tbls-clt-full.cfg](../examples/docking-multiple-ambig/docking-multiple-tbls-clt-full.cfg) example the clustering is performed right after the `rigidbody` module, so as to lump together solutions resulting from the application of different sets of restraints.

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.

## general analysis tasks

In the `analysis` directory, a few examples of how to perform general analysis tasks are provided. These workflows aim at showing how HADDOCK3 can be used outside the scope of docking and refinement to analyse structures and models coming from different sources (e.g. Alphafold2). These example workflows include:

- [alascan-test.cfg](../examples/analysis/alascan-test.cfg): an example of how to perform alanine scanning on a protein-protein complex.
- [contmap-test.cfg](../examples/analysis/contmap-test.cfg): an example of how to generate contact maps from a set of models.
- [topoaa-caprieval-test.cfg](../examples/analysis/topoaa-caprieval-test.cfg): an example of how to run a CAPRI-like (`caprieval`) analysis on an ensemble of models.
- [topoaa-clustfcc-test.cfg](../examples/analysis/topoaa-clustfcc-test.cfg): an example of how to run a FCC clustering analysis on an ensemble of models.
- [topoaa-ilrmsdmatrix-clustrmsd-test.cfg](../examples/analysis/topoaa-ilrmsdmatrix-clustrmsd-test.cfg): an example of how to run a ilRMSD matrix clustering analysis on an ensemble of models.
