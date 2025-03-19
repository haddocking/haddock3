# Installation

Create a virtual environment with Python 3.9:

You can use Python's `venv` or `conda` depending on your choice.
Commands are provided below:

### with `venv`

```bash
python -m venv .venv
source .venv/bin/activate
```

### with `conda`

```bash
conda create -n haddock3 python=3.9
conda activate haddock3
```

```bash
git clone https://github.com/haddocking/haddock3.git
cd haddock3
pip install .
```

## (Optional) Install MPI libraries if you intend to run HADDOCK3 with MPI

To use the mpi implementation of haddock3 you must have mpi4py installed in the `(haddock3)` python environment, and OpenMPI in the host system.

```bash
$ pip install mpi4py
# or
$ conda install -c conda-forge mpi4py
```

Later, you can find [here](https://www.bonvinlab.org/haddock3/tutorials/mpi.html) instructions on how to run HADDOCK3 with MPI.

## Troubleshooting the CNS executable

Depending on your architecture, the CNS executable coming with the `pip install` command (see above), might give errors because of some missing libraries.
In such a case you should recompile CNS (see the [`CNS.md`](CNS.md) file in the `docs` directory). Once you have tested your newly compiled executable you can place it in the install haddock3 version. For this do the following:

1. First get the location where haddock3 has been installed, e.g.:

```bash
> pip show haddock3
Name: haddock3
Version: 2024.10.0b7
Summary: HADDOCK3
Home-page:
Author:
Author-email: BonvinLab <bonvinlab.support@uu.nl>
License: Apache License 2.0
Location: /home/testuser/.venv/lib/python3.9/site-packages
...
```

2. Copy your CNS executable to the haddock3 installation:

```bash
> cp <my-cns-binary> /home/testuser/.venv/lib/python3.9/site-packages/haddock/bin/cns
```

## (Optional) Install web service dependencies if you intend to run HADDOCK3 restraints web service

To run the restraints web service you must have the following dependencies installed in the `(haddock3)` python environment:

```bash
pip install uvicorn fastapi
```

Information on the restraints web service can be found [here](https://github.com/haddocking/haddock3/blob/main/src/haddock/clis/restraints/webservice.py).

## Installing third-party packages

HADDOCK3 can integrate third-party software in its workflows.
We are not responsible for the proper installation of such packages, but
we help you install them. Below, you will find a list of all third-party
packages HADDOCK3 can use and guidelines for their proper installation.

### `lightdock` (outdated)

To install [lightdock](https://github.com/lightdock/lightdock) follow
the instructions on the project's website. Remember to install it under
the same Python environment you created for HADDOCK3. If you have any
doubts, please let us know.

### `gdock` (outdated)

1. Clone the latest version:

```
cd some-folder
git clone https://github.com/rvhonorato/gdock.git
```

2. Install Python3+ dependencies

```
pip install deap scipy mgzip biopython
```

3. Set `GDOCK_PATH`

```
export GDOCK_PATH=some-folder
```

**Important**: These are not the full `gdock` installation
instructions as only the model generation is used here. Please check the
[repository page](https://github.com/rvhonorato/gdock) for more
information.

## `openmm`

The use of the openmm module requires the installation of both OpenMM and pdbfixer.

### In venv

The venv installation instruction uses `pip` and requires a python version >=3.10 (not functional with python3.9).

```bash
# Create a new virtual env
python3.12 -m venv hd3-openmm
source hd3-openmm/bin/activate
# Install haddock3
pip install .
# Install OpenMM
pip install openmm==8.2.0
# Install openmm - pdbfixer
pip install https://github.com/openmm/pdbfixer/archive/refs/tags/v1.10.tar.gz
```

### In conda

The Conda installation instructions are using conda-forge to retrieve packages and requires a python3 version between 3.9 and 3.13.

```bash
conda create -n hd3-p39-openmm python=3.9
# Activate the haddock3 env
conda activate hd3-p39-openmm
# Install haddock3
pip install haddock3
# libstdxx must be installed, maybe already present in your system
conda install -c conda-forge libstdcxx-ng
# Install OpenMM and related tools/binaries
conda install -c conda-forge openmm==8.2.0 pdbfixer==1.10
```

OpenMM should automatically detect the fastest platform among those available
on your machine.

Please refer to the [official page](http://docs.openmm.org/latest/userguide/)
of the project for a full description of the installation procedure.
