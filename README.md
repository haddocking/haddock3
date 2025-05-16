# `haddock3`

![haddock3-logo](https://raw.githubusercontent.com/haddocking/haddock3/refs/heads/main/docs/figs/HADDOCK3-logo.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10527751.svg)](https://doi.org/10.5281/zenodo.10527751)
[![Research Software Directory](https://img.shields.io/badge/rsd-haddock3-00a3e3.svg)](https://research-software-directory.org/software/haddock3)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/8844/badge)](https://www.bestpractices.dev/projects/8844)

[![Downloads](https://static.pepy.tech/badge/haddock3)](https://pepy.tech/project/haddock3)
![PyPI - Version](https://img.shields.io/pypi/v/haddock3)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/haddock3)

[![ci](https://github.com/haddocking/haddock3/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/haddock3/actions/workflows/ci.yml)
[![pages](https://github.com/haddocking/haddock3/actions/workflows/pages.yml/badge.svg)](https://github.com/haddocking/haddock3/actions/workflows/pages.yml)
[![CodeQL](https://github.com/haddocking/haddock3/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/haddocking/haddock3/actions/workflows/github-code-scanning/codeql)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e11e7f45400f4e8589cdf5941f95233a)](https://app.codacy.com/gh/haddocking/haddock3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e11e7f45400f4e8589cdf5941f95233a)](https://app.codacy.com/gh/haddocking/haddock3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

## Introduction

HADDOCK, standing for **H**igh **A**mbiguity **D**riven protein-protein **DOCK**ing, is a widely used computational tool for the integrative modeling of biomolecular interactions. Developed by researchers at [Utrecht University](https://uu.nl) in the [BonvinLab](https://bonvinlab.org) for more than 20 years, it integrates various types of experimental data, biochemical, biophysical, bioinformatic prediction and knowledge to guide the docking process.

## Installation

Simple installation of the [latest release](https://pypi.org/project/haddock3/) of HADDOCK3 (assuming you have a Python version 3.9 to 3.12 installed and the rights to install the software - if not refer to the instructions in [INSTALL.md](docs/INSTALL.md) for using either `venv` or `conda`) . 

```bash
pip install haddock3
```

In case you rather install the latest unreleased version use instead:

```bash
git clone https://github.com/haddocking/haddock3.git
cd haddock3
pip install .
```

For detailed instructions and installation of third-party software, please check [INSTALL.md](docs/INSTALL.md) 

You might also want to check the following utilities:

- [`haddock-restraints`](https://github.com/haddocking/haddock-restraints): Tool to generate restraints to be used in `haddock3`.
- [`haddock-runner`](https://github.com/haddocking/haddock-runner): Tool to run large scale `haddock3` simulations using multiple input molecules in different scenarios
- [`haddock-tools`](https://github.com/haddocking/haddock-tools): Set of useful utility scripts developed over the years by the BonvinLab group members

## Usage

The most basic usage is:

```bash
haddock3 <configuration-file.toml>
```

For help on haddock3 usage:

```bash
$ haddock3 -h
usage: haddock3 [-h] [--restart RESTART] [--extend-run EXTEND_RUN] [--setup] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [-v] recipe

positional arguments:
  recipe                The input recipe file path

optional arguments:
  -h, --help            show this help message and exit
  --restart RESTART     Restart the run from a given step. Previous folders from the selected step onwards will be deleted.
  --extend-run EXTEND_RUN
                        Start a run from a run directory previously prepared with the `haddock3-copy` CLI. Provide the run directory created with `haddock3-copy` CLI.
  --setup               Only setup the run, do not execute
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
  -v, --version         show version
```

Check the [EXAMPLES](https://github.com/haddocking/haddock3/blob/main/examples/README.md) page for more some usage examples and the [MANUAL]() (_coming soon!_) for a more detailed explanation of the configuration file.

## Support

If you encounter any code-related issues, [please open an issue](https://github.com/haddocking/haddock3/issues/new/choose).

If you have any other questions or need help, please contact us at [ask.bioexcel.eu](https://ask.bioexcel.eu/).

If you clone this repository and use `haddock3` for your research, please support us by signing up in [this form](https://forms.gle/LCUHiYHh1hE9rd8L6). This will allow us contact you when needed for `haddock3`-related issues, and also provide us a mean to demonstrate impact when reporting for grants - which grealty helps us to keep the project alive!

## Cite us

If you used `haddock3` for your research, please cite us:

- **Research article**: M. Giulini, V. Reys, J.M.C. Teixeira, B. Jiménez-García, R.V. Honorato, A. Kravchenko, X. Xu, R. Versini, A. Engel, S. Verhoeven, A.M.J.J. Bonvin, *HADDOCK3: A modular and versatile platform for integrative modelling of biomolecular complexes* biorxiv (2025). [https://www.biorxiv.org/content/10.1101/2025.04.30.651432v1](https://www.biorxiv.org/content/10.1101/2025.04.30.651432v1)
 
- **Cite this repository**: M.C. Teixeira, J., Vargas Honorato, R., Giulini, M., Bonvin, A., SarahAlidoost, Reys, V., Jimenez, B., Schulte, D., van Noort, C., Verhoeven, S., Vreede, B., SSchott, & Tsai, R. (2024). haddocking/haddock3: v3.0.0-beta.5 (Version 3.0.0-beta.5) [Computer software]. [https://doi.org/10.5281/zenodo.10527751](https://doi.org/10.5281/zenodo.10527751)

## Useful resources

- [User manual]() (_coming soon!_)
- [Best practice guide](https://www.bonvinlab.org/software/bpg/)
- [The HADDOCK2.4 web server: A leap forward in integrative modelling of biomolecular complexes. Nature Prot. 2024](https://www.nature.com/articles/s41596-024-01011-0)

## Development

Please check [DEVELOPMENT](https://github.com/haddocking/haddock3/blob/main/DEVELOPMENT.md) for instructions on how to develop `haddock3`

### Code Documentation

The code documentation is automatically built and hosted at [bonvinlab.org/haddock3](https://www.bonvinlab.org/haddock3/).

To build it locally it is necessary to have some extra packages installed. You can install them using the following command:

```bash
pip install -e '.[docs]'
```

Then, to build the documentation, run the following commands:

```bash
sphinx-apidoc -f -e -o docs/ src/haddock -d 1
sphinx-build -b html docs haddock3-docs
```

> Warning messages are expected, but the documentation should be built successfully.

The rendered documentation will be available at `haddock3-docs/index.html`. This will open a local webpage with the

### Contributing

Check the [CONTRIBUTING](CONTRIBUTING.md) file for instructions on how to contribute with the project!

<!-- ---

Happy HADDOCking!

<img src="https://www.bonvinlab.org/images/bio-haddock.png" alt="haddock" width="50px"> -->
