[![Build Status](http://alembick.science.uu.nl:8080/buildStatus/icon?job=haddock3%2Fmaster)](http://alembick.science.uu.nl:8080/job/haddock3/job/master/)
[![codecov](https://codecov.io/gh/haddocking/haddock3/branch/master/graph/badge.svg?token=K2UshyxoRu)](https://codecov.io/gh/haddocking/haddock3)

![HADDOCK3](docs/media/HADDOCK3-logo.png)


**ATTENTION: This repository is under heavy development and may change abruptly.**

# Changelog

* 0.0.alpha1 (06-11-2019)
    * First version of the skeleton code, 
    * Simple protein-protein with ambiguious restraints
    * Scoring of an ensamble of models (clustering included)

# TODO

**WIP**

# Instalation

* Requirements
    * [CNS 1.31 UU](https://surfdrive.surf.nl/files/index.php/apps/files/?dir=/Shared/HADDOCK/CNS&fileid=5041663829)
    * Python 3.7.x

```bash
$ git clone https://github.com/haddocking/haddock3.git
$ cd haddock3
$ setenv PYTHONPATH ${PYTHONPATH}:`pwd`
$ pip install -r requirements.txt --user
$ cd haddock/src
$ make
$ chmod +x contact_fcc
$ cd ../../

# Edit "cns_exe" and "haddock3" in the ini script
$ vim haddock/etc/haddock3.ini
```

# Execution

```bash
$ cd examples/protein-protein
$ python3 ../../haddock/haddock3.py run.toml
```

# Scoring

```bash
$ cd examples/protein-protein
$ python ../../haddock/workflows/scoring/setup_scoring.py scoring.toml
$ cd run-scoring-example
$ python3 ../../../haddock/workflows/scoring/run_scoring.py
```


***
# Development

The default recipes are refactored versions of the "legacy" 
protocols (`generate.inp`, `refine.inp`, `refine_h2o.inp` and `scoring.inp`) now called: 

   * `generate-topology.cns`
   * `it0.cns`
   * `it1.cns`
   * `itw.cns`
   * `scoring.cns`

These recipes are independent from each other and compatible with the "modular" manner they are executed, 
example: 

By using `@RUN:protocols/scale_inter_final.cns`, the `scale_inter_final.cns` script is executed alongside the main 
`CNS` instance, having access to all variable definitions. To emulate this we developed a Recipe module 
(`haddock.modules.worker.recipe`) that reads the main `CNS` code, identify its dependency tree and recursively append
 each of the scripts to the appropriate position. 
 
 This step is done by `setup_haddock.py`, which also adjusts the user-specified parameters, saving a single template 
  to, for example, `run1/topology/template/` resulting in a longer, but self-contained, `.inp`. With this design 
  there is no need to copy over `protocols/` and `toppar/` to the simulation folder, reducing I/O and taking less space.

Having the execution separated in two steps allows the user to manually check the simulation or to prepare a large 
batch of files. However manually editing the templates, as one would edit `run.cns` is possible but should no longer be 
necessary. Each recipe has a companion `.json` that defines its default parameters (`it0.cns`/`it0.json`), 
non-default parameters are passed by the user via the `run.toml` file, an upgrade of `run.param` (or `new.html`).

 TOML ([Tom's Obvious, Minimal Language.](https://github.com/toml-lang/toml])) was chosen since it is human readable 
 (allowing for comments) and an efficient low-impact way of passing parameters, example `run.toml`:
 
```toml
title = "HADDOCK3 Setup file"
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
sampling = 200
params.auto_his = true
params.noecv = false

[stage.semi_flexible]
recipe='default'
sampling = 20

[stage.water_refinement]
recipe='default'
sampling = 20
#===========================================================#
```

In this example the simulation scheme (reminiscent of the `Queue`) will be parallel, `it0` will use its default 
parameters with the exception of `auto_his` and `noecv` which were manually defined with `sampling=200`, `it1` and 
`itw` will be default with `n=20`.

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


**WIP**

## CNS Refactoring guidelines

### From global to local variables


The `noecv` parameter, which is responsible for random removal of restraints. In `run.cns`, `noecv` is 
defined and evaluated as follows, becoming a global variable

```
define (...{===>} noecv=true;...)
...
evaluate (&data.noecv=&noecv)
```

This parameter is then called in `protocols/refine.inp`, first by reading it and then making it usable in the rest of 
the code as 
`$Data.noecv`.

```
@RUN:run.cns(...Data=$Data;...)
...
if ($Data.noecv eq true) then ...
```

To get the same in the new implementation, the `CNS` code **must** be refactored. Variables are now defined 
"*on-the-fly*" by a combination of reading a `.json` with default values and a `.toml` with user custom parameters 
instead directly editing `run1/run.cns`

```json
{
  "params": {
    "noecv": true
  }
}
```

The issue arrives because "legacy" `CNS` protocols will use `$Data.noecv` instead of `$noecv`, to bypass this the 
main `CNS` protocol of a given recipe **needs** to be manually edited:
```
evaluate ($Data.noecv=$noecv)
```

New `CNS` protocols should rely on the literal name of the variable, `$noecv` instead of `$Data.noecv`.
