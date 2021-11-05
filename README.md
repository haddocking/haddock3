# WARNING: The `main` branch is not production-ready

The `main` branch is a prototype of newly refined architecture and it
does not yet contain the functionalities we have reported previously.
For a running version of HADDOCK3 please refer to the `alpha1` or the
release page. However, we won't develop `alpha1` further. Stay tuned for
new updates on the `main` branch as we are actively working on it.
Cheers!

* * *

# HADDOCK3
## 1. Installation

### 1.1 Clone this repository:

```bash
git clone --recursive https://github.com/haddocking/haddock3.git
cd haddock3
cd src/fcc/src
chmod u+x Makefile
./Makefile
cd -
```

### 1.2 Create a virtual environment with Python 3.8+ and install dependencies:
#### with `venv`

```bash
virtualenv-3.8 venv
source venv/bin/activate
pip install -r requirements.txt
python setup.py develop --no-deps
```

#### with `conda`
```bash
conda env create -f requirements.yml
conda activate haddock3
python setup.py develop --no-deps
```

### 1.3 Make a CNS binary shortcut to the expected path:

```bash
mkdir -p bin/

# on mac
ln -s /PATH/TO/cns_solve-1.31-UU-MacIntel.exe bin/cns

# on linux
ln -s /PATH/TO/CNS_FOLDER/intel-x86_64bit-linux/source/cns_solve-2002171359.exe bin/cns
```

### 1.4 Keep up to date

In the `github` folder of `haddock3` run:

```bash
git pull --recurse-submodules
```

This will pull the latest changes to your local folder and because you
installed `haddock3` with the option `develop` those changes become
available immediately.

## 2. Examples

### 2.1. Basic scoring of an ensemble of 5 structures:

```bash
cd examples/recipes/scoring/
haddock3 scoring.toml
```
