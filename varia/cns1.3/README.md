### CNS-related code and scripts

HADDOCK requires that CNS1.3 be recompiled to increase in particular the size of some arrays and adding new functionalities.
For this, when installing CNS, copy first all the files from this directory into the source directory of your CNS installation before starting the installation.
The large number of molecule types and their corresponding force field parameters does require to increase the size of various arrays.

To copy the code and scripts to your CNS installation directory use:

```bash
cp -r source bin instlib <CNS-installation-directory>
```

The provided getarch script under the bin dir adds support for Mac M1/2/3 processors.

The provided Makefiles under the instlib dir adds the proper gfortran compiler options.

And the provided code adds functionality for:


#### Radius of gyration restraints

A radius of gyration distance restraint can be used in HADDOCK. By default it is applied to the entire system, but can be restricted to part of the system using standard CNS atom selections. For details refer to the online [HADDOCK manual](https://www.bonvinlab.org/software/haddock2.4/Rg/).

_Note_ that this functionality is not yet implemented in HADDOCK3


#### Residual dipolar couplings defined as intervector projection angle restraints

Intervector projection angle restraints ( [Meiler et al. _J. Biomol. NMR_ **17**, 185 (2000)](https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=10805131&dopt=Abstract)) are obtained by taking pairs of residual dipolar couplings and generating intervector projection angle restraints (somewhat equivalent to dihedral angle restraints). These restraints have the advantage that they do no longer depend on the orientation of the dipole vector with respect to the alignment tensor. Instead they restrain the angle between two dipolar vectors, allowing for two minima. Two force constants must be therefore defined: one for the border potential function and one for the central part (e.g. between the two minima). For details refer to the online [HADDOCK manual](https://www.bonvinlab.org/software/haddock2.4/RDC/#intervector-projection-angle-restraints-for-docking).


_Note_ that this functionality is not yet implemented in HADDOCK3


#### Pseudo contact shift restraints

Pseudo contact shifts (PCS) can provide useful information on both the distances and the orientation of the molecules to be docked. They can be used directly as restraints in HADDOCK using the XPCS energy term in CNS. For details refer to the online [HADDOCK manual](https://www.bonvinlab.org/software/haddock2.4/PCS/).


_Note_ that this functionality is not yet implemented in HADDOCK3

