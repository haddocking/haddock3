[![Build Status](http://alembick.science.uu.nl:8080/buildStatus/icon?job=haddock3%2Fmaster)](http://alembick.science.uu.nl:8080/job/haddock3/job/master/)
[![codecov](https://codecov.io/gh/haddocking/haddock3/branch/master/graph/badge.svg?token=K2UshyxoRu)](https://codecov.io/gh/haddocking/haddock3)

![HADDOCK3](docs/media/HADDOCK3-logo.png)

The official repo of the new modular BioExcel2 version of HADDOCK.

**ATTENTION: This repository is under heavy development and may change abruptly.**

***
## Stages

1. Initialization
2. Pre-processing
3. Topology generation
4. Docking
    1. Rigid-body
    2. Semi-Flexible
    3. Water-refinement

### 1. Initialization
### 2. Pre-processing
### 3. Topology generation
### 4. Docking

## Workflows

* Scoring
* Refinement

***
# Dev information

## CNS Refactoring guidelines

### Variables
The stable implementation of HADDOCK relies heavily on global variables.

Example: `noecv`  

In `run.cns` the `noecv` variable is defined and evaluated to its global variable:

```
define (
...
{===>} noecv=true;
...
)
...
evaluate (&data.noecv=&noecv)
...
```

We can then see in `refine.inp` how this variable is used. First by reading it and make it usable as `$Data.noecv` anywhere in the code.

```
@RUN:run.cns(
...
Data      =$Data;
...
)
...
if ($Data.noecv eq true) then
...
```

To achieve this behavious in a modular sense, part of the CNS code must be refactored to account for the new input 
method. In the new 
implementation variables will be defined on the fly, reading from a `json` parameter file instead of `run.cns` and 
written on the header of the input file.

Example:
```
# Json file
{
  "params": {
    "noecv": true
  }
}

# CNS input file header
evaluate ($noecv=true)
```

However `refine.inp` and its dependecies use `$Data.noecv` instead of `$noecv`, so this must be manually accounted 
for simply by adding the following to the top of the refactored CNS recipe.
```
evaluate ($Data.noecv=$noecv)
```

Future edits should not use the `$Data.` format and stick to the literal variable instead.
***
