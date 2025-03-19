# Developing HADDOCK3

This file provides information on how to setup a development environment for HADDOCK3.

- [System requirements](#system-requirements)
- [Setting up the development environment](#setting-up-the-development-environment)
  - [Clone the repository](#clone-the-repository)
  - [Python environment](#python-environment)
  - [Install haddock3 in development mode](#install-haddock3-in-development-mode)
  - [Running tests](#running-tests)

## System requirements

- Python 3.9
- OpenMPI

### Installing system dependencies

Below the instructions are provided for a Ubuntu system using `apt-get`. If you are using a different system, please refer to the respective package manager - `yum`, `dnf`, `pacman`, `homebrew`, `port`, etc.

```bash
sudo apt-get update &&
    # Needed for building Python
    sudo apt-get install build-essential python-setuptools python-pip &&
    sudo apt-get install libncursesw5-dev libgdbm-dev libc6-dev &&
    sudo apt-get install zlib1g-dev libsqlite3-dev tk-dev &&
    sudo apt-get install libssl-dev openssl &&
    sudo apt-get install libffi-dev &&
    # Needed for haddock3 installation
    sudo apt-get install git gcc &&
    # Needed to run the MPI related tests
    sudo apt-get install openmpi-bin libopenmpi3 libopenmpi-dev
```

After all the system-dependencies are in place, download and compile python.

**Conda is NOT recommended for development**

```bash
wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz
tar -xf Python-3.9.6.tgz
cd Python-3.9.6
./configure --enable-optimizations
sudo make altinstall -j 8
```

Then `python3.9` should be available on your system at `/usr/local/bin/python3.9`

On **OSX**, you can use a package-manager such as [brew](https://brew.sh) to install Python 3.9.

> Please keep in mind installing python with a package manager can mask system dependencies that your development might add. That's why we recommended you install it from source.

```bash
brew install python@3.9
```

## Setting up the development environment

### Clone the repository

```bash
git clone https://github.com/haddocking/haddock3.git
cd haddock3
```

### Python environment

We recommend you use Python's native virtual environment to manage the dependencies.

```bash
python3.9 -m venv .venv
source .venv/bin/activate
```

### Install haddock3 in development mode

Install both project dependencies and test dependencies using pip.

```bash
pip install -e '.[dev,docs]'
```

> If you are using a Mac, if the installation of mpi4py fails, run first `brew install mpi4py`

## Running tests

In `haddock3` we use the pytest framework, make sure you check it's [documentation](https://docs.pytest.org/en/6.2.x/contents.html) the tests are located in `tests/` (unit) and `integration_tests/` directories.

```bash
pytest tests/
```

```bash
pytest integration_tests/
```

## Troubleshooting the CNS executable

Depending on your architecture, the CNS executable coming with the `pip install` command (see above), might give errors because of some missing libraries.
In such a case you should recompile CNS (see the [`CNS.md`](docs/CNS.md) file in the `docs` directory). Once you have tested your newly compiled executable you can place it in the install haddock3 version. For this do the following:

1) First get the location where haddock3 has been installed:

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

2) Copy your CNS executable to the haddock3 installation:

```bash
> cp <cns.binary> /home/testuser/.venv/lib/python3.9/site-packages/haddock/bin/cns
```

## Installation in an HPC environment

**Please get in contact with the system administrator before doing development in a shared HPC environment.**

For installation in an HPC environment we recommend to check the installed Python versions on the system and also importantly if an `openmpi` (or other custom MPI) installation is available on the system.
Those are often offered via the `module` command.
If you only intend to develop haddock3 using the multiprocessing scheduler, the above instructions should be fine. But to harvest the MPI capabilities of an HPC system it is best to build haddock3 using the installed MPI version on the HPC system.
