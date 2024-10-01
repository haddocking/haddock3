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

1. Install the latest version of openmm and pdbfixer
```bash
# Activate the haddock3 env
conda activate haddock3
# libstdxx must be installed, maybe already present in your system
conda install -c conda-forge libstdcxx-ng
# Install OpenMM and related tools/binaries
conda install -c conda-forge openmm pdbfixer
```

Openmm should automatically detect the fastest platform among those available
on your machine.

Please refer to the [official page](http://docs.openmm.org/latest/userguide/)
of the project for a full description of the installation procedure.
