# Installation

Open a `terminal` window and navigate to the folder where you want to
install HADDOCK3; for example: `software`. The current installation
instructions are local and will affect only your user.

Before starting with the installation of HADDOCK3, make sure to properly [install CNS](CNS.md).
If you have installed a previous version of HADDOCK, you may already have a suitable version of CNS.
Please do [check your CNS installation](CNS.md#5-Check-installation) before proceeding.


## 1. Clone this repository:

Mind the `--recursive` flag when cloning!

```bash
git clone --recursive https://github.com/haddocking/haddock3.git
cd haddock3
cd src/fcc/src
chmod u+x Makefile
make
cd -
```

By the end of the above commands, you should be back to the `haddock3`
main folder.

## 2 Create a virtual environment with Python 3.9+ and install dependencies:

You can use Python's `venv` or `conda` depending on your choice.
Commands are provided below:

### with `venv`

```bash
virtualenv venv --python=3.9
source venv/bin/activate
pip install -r requirements.txt
```

### with `conda`

```bash
conda create -n haddock3 python=3.9
conda activate haddock3
pip install -r requirements.txt
```

## 3. Install the HADDOCK3 package and command line clients

```bash
python setup.py develop --no-deps
```

## 4. Make a CNS binary shortcut to the expected path:

```bash
mkdir -p bin/

# on mac
ln -s /PATH/TO/cns_solve_1.3/mac-intel-darwin/source/cns_solve-2206031450.exe bin/cns

# on linux
ln -s /PATH/TO/cns_solve_1.3/intel-x86_64bit-linux/source/cns_solve-2002171359.exe bin/cns
```

As long as you have the HADDOCK3 python environment activated, you can
navigate away from the HADDOCK3 installation folder. You can run
HADDOCK3 from anywhere. To run HADDOCK3, follow the [usage
guidelines](USAGE.md).


## 5. Keep your installation up to date

Navigate to the `haddock3` installation folder (the one you cloned from
GitHub). Ensure you have the `haddock3` python environment activated.
Please keep in mind that HADDOCK3, as well as its
dependencies, are under active development.
If the updating processing fails, it is safest to reinstall
from scratch. Always refer to the latest installation guidelines.

```bash
# if you used `venv`
source venv/bin/activate

# if you used `conda`
conda activate haddock3
```

Afterwards:

```bash
# pull the latest source code from our repository to your computer
git pull

# if you used venv to create the python environment, run:
pip install -r requirements.txt  --upgrade

# if you used anaconda to create the python environment, run:
conda env update -f requirements.yml

# ensure all command-lines clients are installed
python setup.py develop --no-deps
```


## 6. (Optional) Install MPI libraries if you intend to run HADDOCK3 with MPI

To use the mpi implementation of haddock3 you must have mpi4py installed in the haddock3 python environment, and OpenMPI in the host system.

```bash
$ pip install mpi4py
# or
$ conda install -c conda-forge mpi4py
```

Later, you can find [here](https://www.bonvinlab.org/haddock3/tutorials/mpi.html) instructions on how to run HADDOCK3 with MPI.

# Installing third-party packages

HADDOCK3 can integrate third-party software in its workflows.
We are not responsible for the proper installation of such packages, but
we help you install them. Below, you will find a list of all third-party
packages HADDOCK3 can use and guidelines for their proper installation.

## `lightdock`

To install [lightdock](https://github.com/lightdock/lightdock) follow
the instructions on the project's website. Remember to install it under
the same Python environment you created for HADDOCK3. If you have any
doubts, please let us know.

## `gdock`

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
