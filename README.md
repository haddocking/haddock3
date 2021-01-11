# HADDOCK3

## 1. Installation

* Clone this repository:

```bash
git clone https://github.com/brianjimenez/haddock3.git
cd haddock3
```

* Create a virtual environment with Python 3.8+ and install dependencies:

```bash
virtualenv-3.8 venv
source venv/bin/activate
pip install -r requirements.txt
```

* Activate and copy CNS binary to the expected path:

```bash
source bin/activate_haddock
mkdir -p bin/cns
cp /path/to/cns_solve-1.31-UU-MacIntel.exe bin/cns
```

* Define some addiational environmental variables:

```bash
export HADDOCK3_NUM_CORES=8
export HADDOCK3_CNS_EXE="/path/to/haddock3/bin/cns/cns_solve-1.31-UU-MacIntel.exe"
```


## 2. Examples

### 2.1. Basic scoring of an ensemble of 5 structures:

```bash
cd examples/recipes/scoring/
haddock3.py scoring.toml
```

