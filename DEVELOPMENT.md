# Developing HADDOCK3

This file provides information on how to setup a development environment for HADDOCK3.

- [System requirements](#system-requirements)
  - [Python3.9](#python39)
  - [CNS](#cns)
  - [git](#git)
  - [gcc](#gcc)
  - [OpenMPI](#openmpi)
- [Setting up the development environment](#setting-up-the-development-environment)
  - [Clone the repository](#clone-the-repository)
  - [Link CNS binary](#link-cns-binary)
  - [Install shipped C++ dependencies](#install-shipped-c-dependencies)
    - [fcc](#fcc)
    - [fast-rmsdmatrix](#fast-rmsdmatrix)
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

Conda is NOT recomended for development. Instead, install and compile Python 3.9 from source.

```bash
wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz
tar -xf Python-3.9.6.tgz
cd Python-3.9.6
./configure --enable-optimizations
sudo make altinstall -j 8
```

Then `python3.9` should be available on your system, at `/usr/local/bin/python3.9`

#### CNS

Please refer to the [CNS installation instructions](docs/CNS.md).

#### Git

```bash
sudo apt-get install git
```

#### gcc

```bash
sudo apt-get install gcc
```

#### OpenMPI

```bash
sudo apt-get install openmpi-bin libopenmpi3 libopenmpi-dev
```

## Setting up the development environment

### Clone the repository

```bash
git clone --recursive https://github.com/haddocking/haddock3.git
```

### Link CNS binary

```bash
mkdir haddock3/bin
ln -s /path/to/cns_solve_1.3/bin/cns haddock3/bin/cns
```

### Install shipped C++ dependencies

#### fcc

```bash
cd src/fcc/src
make
```

#### fast-rmsdmatrix

```bash
cd src/fast-rmsdmatrix/src
make fast-rmsdmatrix
```

### Python environment

We recomend you use Python's native virtual environment to manage the dependencies.

```bash
/usr/local/bin/python3.9 -m venv .venv
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

## Running tests

Simply run the following command to run the tests.

```bash
pytest
```
