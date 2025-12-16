# Installation

Create a virtual environment with Python3.9+

You can use Python's `venv` or `conda` depending on your choice.
Commands are provided below:

## with `venv`

```bash
python -m venv .venv
source .venv/bin/activate
```

## with `conda`

```bash
conda create -n haddock3 python=3.9
conda activate haddock3
```

```bash
pip install haddock3
```

## (Optional) Install MPI dependencies

```bash
pip install 'haddock3[mpi]'
```

> Instruction on how to use MPI with haddock3 are found in the [tutorials](https://www.bonvinlab.org/haddock3/tutorials/mpi.html)

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

## `voroscoring`

To be able to use the `[voroscoring]` module, several steps are required:
- haddock3 **MUST** be installed under a conda env
- A functional installation of [FTDMP](https://github.com/kliment-olechnovic/ftdmp)
- The setup of a conda environement (e.g.: ftdmp), in which you will install FTDMP (see: [Help in installing a functional installation of ftdmp](#help-in-installing-a-functional-installation-of-ftdmp))

Once those three conditions are fulfilled, when using the `[voroscoring]` module in haddock3, the configuration file must be tuned to contain parameters describing how to load the appropriate conda env (ftdmp) and where to find FTDMP scripts and executables:

```toml
[voroscoring]
# This parameter defines the base directory where conda / miniconda is installed
conda_install_dir = "/absolute/path/to/conda/"
# This parameter defines the name of the conda env that you created and where FTDMP is installled
conda_env_name = "ftdmp"
# This parameter defines where FTDMP scripts / executables can be found
ftdmp_install_dir = "/absolute/path/to/FTDMP/"
```

### Help in installing a functional installation of ftdmp

Installing ftdmp under a conda environnement can be difficult.

Here are the steps we performed to install it on our cluster:

```bash
conda create -n ftdmp python=3.11
conda activate ftdmp
pip install torch==2.5.0 torchvision==0.20.0 torchaudio==2.5.0 --index-url https://download.pytorch.org/whl/cu124
python -c "import torch; print(torch.__version__)"
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.5.0+cu124.html
pip install torch_geometric
conda install -c conda-forge openmm pdbfixer libstdcxx-ng
conda install r r-essentials --channel conda-forge
```
