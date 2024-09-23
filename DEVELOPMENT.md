# Developing HADDOCK3

This file provides information on how to setup a development environment for HADDOCK3.

- [System requirements](#system-requirements)
  - [Python3.9](#python39)
  - [git](#git)
  - [gcc](#gcc)
  - [OpenMPI](#openmpi)
- [Setting up the development environment](#setting-up-the-development-environment)
  - [Clone the repository](#clone-the-repository)
  - [Python environment](#python-environment)
  - [Install Python dependencies](#install-python-dependencies)
  - [Running tests](#running-tests)

## System requirements

- Python 3.9
- CNS
- git
- gcc
- OpenMPI

### Installing system dependencies

Below the instructions are provided for a Ubuntu system using `apt-get`. If you are using a different system, please refer to the respective package manager - `yum`, `dnf`, `pacman`, `homebrew`, `port`, etc.

#### Python3.9

**Conda is NOT recommended for development. Instead, install and compile Python 3.9 from source.**

On a **Linux-based system** you can install Python 3.9 with the following commands:

```bash
wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz
tar -xf Python-3.9.6.tgz
cd Python-3.9.6
./configure --enable-optimizations
sudo make altinstall -j 8
```

Then `python3.9` should be available on your system at `/usr/local/bin/python3.9`

On **OSX**, use rather [brew](https://brew.sh) to install Python 3.9:

```bash
brew install python@3.9
```

Then `python3.9` should be available on your system at `/opt/homebrew/bin/python3.9`

#### Git

On a **Linux-based** system:

```bash
sudo apt-get install git
```

#### gcc

On a **Linux-based** system:

```bash
sudo apt-get install gcc
```

#### OpenMPI

On a **Linux-based** system:

```bash
sudo apt-get install openmpi-bin libopenmpi3 libopenmpi-dev
```

**Note** that this step is not required on OSX.

## Setting up the development environment

### Clone the repository

```bash
git clone --recursive https://github.com/haddocking/haddock3.git
```

### Python environment

We recomend you use Python's native virtual environment to manage the dependencies.

```bash
/usr/local/bin/python3.9 -m venv .venv
source .venv/bin/activate
```

On **OSX** the command would be:

```bash
/opt/homebrew/bin/python3.9 -m venv .venv
source .venv/bin/activate
```

### Install Python dependencies

Install both project dependencies and test dependencies using pip.

```bash
pip install -r requirements.txt &&
  pip install \
    coverage==7.2.5 \
    pytest==7.3.1 \
    pytest-cov==4.0.0 \
    hypothesis==6.75.1 \
    pytest-mock==3.12.0 \
    fastapi==0.110.1 \
    httpx==0.27.0 \
    mpi4py==3.1.6
```

Install haddock3 in development mode.

```bash
python setup.py develop
```

## Installation in an HPC environment

**Please get in contact with the system administrator before doing development in a shared HPC environment.**

For installation in an HPC environment we recommend to check the installed Python versions on the system and also importantly if an `openmpi` (or other custom MPI) installation is available on the system.
Those are often offered via the `module` command.
If you only intend to develop haddock3 using the multiprocessing scheduler, the above instructions should be fine. But to harvest the MPI capabilities of an HPC system it is best to build haddock3 using the installed MPI version on the HPC system.

## Running tests

In `haddock3` we use the pytest framework, the tests are located in `tests/` (unit) and `integration_tests/` directories.

```bash
pytest
```
