# `haddock3 notebooks`

![haddock3-logo](https://raw.githubusercontent.com/haddocking/haddock3/refs/heads/main/docs/figs/HADDOCK3-logo.png)

This directory contains Jupyter notebooks that can be directly launched on Google Colab.

To run it locally on your system see the instructions below (does require a working python3 (3.9 to 3.13) installation).

---

## Jupyter notebooks

| Notebooks                           |      Description                     |    Colab    |
|-------------------------------------|--------------------------------------|--------------|
| HADDOCK3-antibody-antigen.ipynb     | antibogy-antigen tutorial (based on our [online tutorial](https://www.bonvinlab.org/education/HADDOCK3/HADDOCK3-antibody-antigen/)) | [Launch Colab](https://colab.research.google.com/github/haddocking/haddock3/blob/main/notebooks/HADDOCK3-antibody-antigen.ipynb) |

---

## Instructions for local execution

```bash
# create a directory
cd $HOME
mkdir haddock3-tutorial
cd haddock3-tutorial

# setup python
python3.13 -m venv .venv
source .venv/bin/activate

# install jupyter
pip install notebook

# download the notebook
wget https://raw.githubusercontent.com/haddocking/haddock3/refs/heads/main/notebooks/HADDOCK3-antibody-antigen.ipynb

# run the jypyter notebook server
jupyter notebook
```
go to <http://localhost:8888/notebooks/HADDOCK3-antibody-antigen.ipynb> - click `Run All` 

---

## Marimo notebooks

[Marimo](https://marimo.io) notebooks are reactive Python notebooks that run locally. Unlike Jupyter, every cell updates automatically when its inputs change, making them well suited for interactive analysis workflows.

| Notebook | Description |
|----------|-------------|
| `HADDOCK3-interface-analysis.py` | Energy minimization → contact map (chord chart + heatmap) with optional alanine scanning, all driven from a single configuration panel |

### Running the interface analysis notebook

**Prerequisites:** HADDOCK3 installed in development mode (`pip install -e '.[dev]'`) and CNS available on your `PATH` (see `docs/CNS.md`).

```bash
# install marimo (if not already present)
pip install marimo

# launch the notebook from the repository root
marimo run notebooks/HADDOCK3-interface-analysis.py
```

Alternatively you can start it in edit mode to see the Python code and be able to edit it.

```
marimo edit notebooks/HADDOCK3-interface-analysis.py
```

The notebook opens in your browser at `http://localhost:2718`. Upload a multi-chain PDB file, adjust the configuration panel, and click **Run HADDOCK3 Workflow**.

---

## Other Jupyter execution options

### [https://notebooks.egi.eu](https://notebooks.egi.eu)

Free for non-profit, but does require registration with the EGI SSO.

Once logged in and the session is active, start a Jypyter notebook and import the notebook from an URL specifying the same https address as above to the local execution.
