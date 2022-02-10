# 1. Installation

Navigate to the folder where you want to install the HADDOCK3 subfolder.

You can safely follow the commands indicated below for the different
steps.

## 1.1 Clone this repository:

```bash
git clone https://github.com/haddocking/haddock3
cd haddock3/src
git clone https://github.com/haddocking/fcc
cd fcc
git checkout 3a1626de3366f6db8b45caed3bd9fbc6b2881286
cd src
chmod u+x Makefile
make
cd ../../..
```

## 1.2 Create a virtual environment with Python 3.9 and install dependencies:

You can use Python's `venv` or Anaconda depending on your choice.
Commands are provided below:

### with `venv`

```bash
virtualenv venv --python=3.9
source venv/bin/activate
pip install -r requirements.txt

# install the HADDOCK3's Python shell and command-line clients (CLIs) on
# the newly created environment.
python setup.py develop --no-deps
```

### with `conda`
```bash
conda env create -f requirements.yml
conda activate haddock3

# install the HADDOCK3's Python shell and command-line clients (CLIs) on
# the newly created environment.
python setup.py develop --no-deps
```

## 1.3 Make a CNS binary shortcut to the expected path:

```bash
mkdir -p bin/

# on mac
ln -s /PATH/TO/cns_solve-1.31-UU-MacIntel.exe bin/cns

# on linux
ln -s /PATH/TO/CNS_FOLDER/intel-x86_64bit-linux/source/cns_solve-2002171359.exe bin/cns
```

As long as you have the HADDOCK3 python environment activated you can
navigate away from the HADDOCK3 installation folder. You can run
HADDOCK3 from anywhere. To run HADDOCK3, follow the [usage
guidelines](USAGE.md).

## 1.4 Keep up to date

Navigate to the `github` folder of `haddock3`. Ensure you have the
haddock3 python environment activated:

```bash
# if you used `venv`
source venv/bin/activate

# if you used `conda`
conda activate haddock3
```

Afterwards:

```bash
git pull
pip install -r requirements.txt  --upgrade  # for venv
conda env update -f requirements.yml  # for conda
python setup.py develop --no-deps
```

This will pull the latest changes to your local folder and because you
installed `haddock3` with the option `develop` those changes become
available immediately.

* * *

# Installing third-party packages

HADDOCK3 is able to integrate third-party software in its workflows.
However, we are not responsible to the proper installation of such
packages, but we do help you to install them. Bellow, you will find a
list of all third-party packages HADDOCK3 can use and guidelines for
their proper installation.

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
