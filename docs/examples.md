## HADDOCK3 examples


The HADDOCK3 examples directory contains various directories and config files
corresponding to different types of complexes, scenarios and data.

Each directory contains both:

- __test__ runs config files (to test the various workflows, make sure the installation works and any changes to the code/scripts have not broken the machinery). Those are setup to run locally.
- __full__ runs config files with recommended parameter settings. Those runs are set up to be executed in HPC mode using slurm (the `...full.cfg`) files. Examples making use of MPI are also provided in some cases, together with an associated job file that should be submitted to the slurm batch system (`...full-mpi.cfg` and `...full-mpi.job`). Make sure to adapt the full config files to your own system.

The following examples are currently provided:

### docking-antibogy-antigen

An antibody-antigen docking example making use only of the knowledge of the hypervarialbles (HV) loops on the antibody to guide the docking. Two different ways of using the knowledge of the HV loop residues are illustrated:

- __CDR-accessible__: using ambiguous distance restraints (AIRs) between the HV loop residues and the solvent-accesible residues of the antigen. Those are defined in the `ambig.tbl` file provided in the `data` directory
- __ranairCDR__: using random AIRs sampled from the HV loops on the antibody and the solvent-accessible residues on the antigen. The random AIRs are used during the rigidbody docking sampling phase and replaced by contact AIRs during refinement.

Since antibodies consists of two separate chains, a few additional unambiguous restraints are defined between the antibody chains to keep them together during the flexible refinement stage. Those are defined in the `unambig.tbl` file in the `data` directory.

Two different protocols/workflows are illustrated:

1) 10000 rigidbody docking models, selection of top500 and flexible refinement + EM of those
2) 10000 rigidbody docking modles, FCC clustering and selection of max 20 models per cluster followed by flexible refinement and EM (config files with `ctl` in their name).

The `caprieval` module is called at various stages during the workflow to assess the quality of the models with respect to the known reference structure.




