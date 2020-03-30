[![Build Status](http://alembick.science.uu.nl:8080/buildStatus/icon?job=haddock3%2Fmaster)](http://alembick.science.uu.nl:8080/job/haddock3/job/master/)
[![codecov](https://codecov.io/gh/haddocking/haddock3/branch/master/graph/badge.svg?token=K2UshyxoRu)](https://codecov.io/gh/haddocking/haddock3)

![HADDOCK3](docs/media/HADDOCK3-logo.png)


**ATTENTION: This repository is under heavy development and will change abruptly.**

***
    
# Instalation

#### Requirements
 * [Crystallography & NMR System (CNS)](http://cns-online.org/v1.3/)
 * Python 3.7.x

#### CNS
 * Make a request on the [CNS website](http://cns-online.org/v1.3/), download the software and copy the contents of `/cns1.3` to CNS's source folder before compilation.

#### HADDOCK

 * Download the [latest release](https://github.com/haddocking/haddock3/releases) *or* clone the repository

```bash
$ wget https://github.com/haddocking/haddock3/archive/v3.0.alpha2.zip
$ unzip haddock3-3.0.alpha2
$ cd haddock3
$ pip install -r requirements.txt
```
 OR

 * Clone the repository, please not that the repository is under constant change 

```bash
$ git clone https://github.com/haddocking/haddock3.git
$ cd haddock3
$ pip install -r requirements.txt
```

#### Third-party
* Edit `haddock3/haddock/etc/haddock3.ini` to account for the third-party software

# Execution

*Edit `CNS_EXE` in `bin/activate` with the correct path*

```bash
$ cd haddock3
$ source bin/activate
$ cd src/haddock3
$ source bin/activate_haddock
$ haddock3.py
```

Example:
```
$ cd examples/protein-protein
$ haddock3.py run.toml
```

Haddock's main input file is a [.toml](https://github.com/toml-lang/toml) containing the parameters that will be used in the simulation.

```toml
#===========================================================#
title = "HADDOCK3 Example setup file"
#===========================================================#
[molecules]
mol1 = '1AY7_r_u.pdb'
mol1.segid = "A"
mol2 = '1AY7_l_u.pdb'
mol1.segid = "B"

[restraints]
ambig = 'ambig.tbl'

[identifier]
run = 1

[execution_parameters]
scheme = 'parallel'
nproc = 2

# Stage specific parameters
[stage]
[stage.topology]
recipe='default'

[stage.rigid_body]
recipe='default'
sampling = 1000
params.auto_his = true

[stage.semi_flexible]
recipe='default'
sampling = 200

[stage.water_refinement]
recipe='default'
sampling = 200
#===========================================================#
```

For a complete list of paramters, click here <under-construction> 

In this example the simulation will be executed in parallel and use the default rigid-body (**it0**) protocol with automatic histidine protonation, default semi-flexible (**it1**) stage and default water-refinement (**itw**).

Running this setup file will create the following folder structure:

```
run1/
|_ data/
    |_ ambig.tbl
    |_ mol1_1.pdb
    |_ mol2_1.pdb
    |_ run.toml
|_ topology/
   |_template/generate-topology.cns
|_ rigid_body/
   |_template/it0.cns
|_ semi_flexible/
   |_template/it1.cns
|_ water_refinement/
   |_template/itw.cns
```

Here we define `one model = one task`, so during execution each `.inp`, `.out` and resulting files will be created 
inside its own folder:

```
run1/topology/
        |_ generate_0000001.inp
        |_ generate_0000001.out
        |_ generate_0000002.inp
        |_ generate_0000002.out
        |_ mol1_1.pdb
        |_ mol1_1.psf
        |_ mol2_1.pdb
        |_ mol2_1.psf
        |_ template/
            |_ generate_topology.cns
```

# Benchmark

HADDOCK v3.0 alpha2 was benchmarked using the [Protein-Protein Docking Benchmark 5](https://github.com/haddocking/BM5-clean) and compared with the current live version. v2.4

The evaluation was done using the true interface of each complex of (4.9 Å) and is expressed in terms of success rate; the ammount of BM5 targets that have at least one docking solution below the specified threhshold whithin a specified subset of solutions ranked by HADDOCK-score.


![BM5](docs/media/haddock3-0-0-alpha2-BM5.png)
***