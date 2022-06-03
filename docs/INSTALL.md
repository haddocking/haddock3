# Installation

Open a `terminal` window and navigate to the folder where you want to
install HADDOCK3; for example: `software`. The current installation
instructions are local and will affect only your user.

## 0 Install CNS

#### Download CNS

Downloading [CNS (Crystallography and NMR System)](http://cns-online.org/v1.3/) requires a license, which academic users can [request for free](http://cns-online.org/cns_request/). You will be sent a password for the download via e-mail.
Using the link provided in the same e-mail, download CNS to the directory in which you want to install it. Navigate to this directory, then uncompress and extract the archive:

```bash
cd ~/software/
gunzip cns_solve_1.3_all.tar.gz
tar -xvf cns_solve_1.3_all.tar
```

Go to the CNS directory:

```bash
cd cns_solve_1.3
```

If you are using csh/tcsh, edit `cns_solve_env` to point to the correct location (the current directory):

```
# CHANGE THE NEXT LINE TO POINT TO THE LOCATION OF THE CNSsolve DIRECTORY

            setenv CNS_SOLVE PATH/TO/cns_solve_1.3
```

If you are using sh/bash, edit `.cns_solve_env_sh`:

```
# CHANGE THE NEXT LINE TO POINT TO THE LOCATION OF THE CNSsolve DIRECTORY

        CNS_SOLVE=PATH/TO/cns_solve_1.3
```

#### Add HADDOCK routines

HADDOCK requires some adjustments in CNS to function properly. Download the [set of routines to add](https://www.dropbox.com/s/wliubqovuusqdvr/cns.tgz?dl=0). Uncompress the archive and move all of its files into the CNS source directory (`PATH/TO/cns_solve_1.3/source/`).

#### Recompile CNS

Make sure you are inside the main `cns_solve_1.3` directory. Compile CNS with the new routines:

```bash
make install

# For csh/tcsh:
source cns_solve_env
# For sh/bash:
source .cns_solve_env_sh
```

## 1 Clone this repository:

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

You can use Python's `venv` or Anaconda depending on your choice.
Commands are provided below:

### with `venv`

```bash
virtualenv venv --python=3.9
source venv/bin/activate
pip install -r requirements.txt
```

### with `conda`

```bash
conda env create -f requirements.yml
conda activate haddock3
```

## 3 Install the HADDOCK3 package and command line clients

```bash
python setup.py develop --no-deps
```

## 4 Make a CNS binary shortcut to the expected path:

```bash
mkdir -p bin/

# on mac
ln -s /PATH/TO/cns_solve-1.31-UU-MacIntel.exe bin/cns

# on linux
ln -s /PATH/TO/CNS_FOLDER/intel-x86_64bit-linux/source/cns_solve-2002171359.exe bin/cns
```

As long as you have the HADDOCK3 python environment activated, you can
navigate away from the HADDOCK3 installation folder. You can run
HADDOCK3 from anywhere. To run HADDOCK3, follow the [usage
guidelines](USAGE.md).


## 5 Keep your installation up to date

Navigate to the `haddock3` installation folder (the one you cloned from
GitHub). Ensure you have the `haddock3` python environment activated.
Please consider HADDOCK3 is under active development, as well as its
dependencies. If the updating processing fails, it is safe to reinstall
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

* * *

# Installing third-party packages

HADDOCK3 can integrate third-party software in its workflows. However,
we are not responsible for the proper installation of such packages, but
we help you install them. Below, you will find a list of all third-party
packages HADDOCK3 can use and guidelines for their proper installation.

## `lightdock`

To install to [lightdock](https://github.com/lightdock/lightdock) follow
the instructions in the project's website. Remember to install it under
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

**Important**: These are not the full `gdock`'s installation
instructions as here only the model generation is used. Please check the
[repository page](https://github.com/rvhonorato/gdock) for more
information.
