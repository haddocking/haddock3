![haddock3-logo](https://raw.githubusercontent.com/haddocking/haddock3/refs/heads/main/docs/figs/HADDOCK3-logo.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20028293.svg)](https://doi.org/10.5281/zenodo.20028293)
[![Research Software Directory](https://img.shields.io/badge/rsd-haddock3-00a3e3.svg)](https://research-software-directory.org/software/haddock3)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/8844/badge)](https://www.bestpractices.dev/projects/8844)

[![Downloads](https://static.pepy.tech/badge/haddock3)](https://pepy.tech/project/haddock3)
![PyPI - Version](https://img.shields.io/pypi/v/haddock3)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/haddock3)
[![SBGrid Badge](https://img.shields.io/badge/haddock3-blue?label=SBGrid&labelColor=%3D)](https://sbgrid.org/software/titles/haddock3)

[![ci](https://github.com/haddocking/haddock3/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/haddock3/actions/workflows/ci.yml)
[![pages](https://github.com/haddocking/haddock3/actions/workflows/pages.yml/badge.svg)](https://github.com/haddocking/haddock3/actions/workflows/pages.yml)
[![CodeQL](https://github.com/haddocking/haddock3/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/haddocking/haddock3/actions/workflows/github-code-scanning/codeql)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e11e7f45400f4e8589cdf5941f95233a)](https://app.codacy.com/gh/haddocking/haddock3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e11e7f45400f4e8589cdf5941f95233a)](https://app.codacy.com/gh/haddocking/haddock3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

## Introduction

HADDOCK, standing for **H**igh **A**mbiguity **D**riven protein-protein **DOCK**ing, is a widely used computational tool for the integrative modeling of biomolecular interactions. Developed by researchers at [Utrecht University](https://uu.nl) in the [BonvinLab](https://bonvinlab.org) for more than 20 years, it integrates various types of experimental data, biochemical, biophysical, bioinformatic prediction and knowledge to guide the docking process.

## Installation

Simple installation of the [latest release](https://pypi.org/project/haddock3/) of HADDOCK3 (requires Python 3.10+):

```bash
pip install haddock3
```

For Python environment setup and installation instructions, see [docs/PYTHON.md](docs/PYTHON.md).

For development installation or to install the latest unreleased version
please refer to [docs/DEVELOPMENT.md](/docs/DEVELOPMENT.md)

## Usage

The most basic usage is:

```bash
haddock3 <configuration-file.toml>
```

Check the [EXAMPLES](https://github.com/haddocking/haddock3/blob/main/examples/README.md) page for more some usage examples and the [User manual](https://www.bonvinlab.org/haddock3-user-manual) for a more detailed explanation of the configuration file.

## Support

If you encounter any code-related issues, [please open an issue](https://github.com/haddocking/haddock3/issues/new/choose).

If you have any other questions or need help, please contact us at [ask.bioexcel.eu](https://ask.bioexcel.eu/).

If you clone this repository and use `haddock3` for your research, please support us by signing up in [this form](https://forms.gle/LCUHiYHh1hE9rd8L6). This will allow us contact you when needed for `haddock3`-related issues, and also provide us a mean to demonstrate impact when reporting for grants - which grealty helps us to keep the project alive!

## Useful resources

- [`haddock-restraints`](https://github.com/haddocking/haddock-restraints): Tool to generate restraints to be used in `haddock3`.
- [`haddock-runner`](https://github.com/haddocking/haddock-runner): Tool to run large scale `haddock3` simulations using multiple input molecules in different scenarios
- [`haddock-tools`](https://github.com/haddocking/haddock-tools): Set of useful utility scripts developed over the years by the BonvinLab group members
- [User manual](https://www.bonvinlab.org/haddock3-user-manual): The online HADDOCK3 guide describing every aspects of the tool.
- [Best practice guide](https://www.bonvinlab.org/software/bpg/) (HADDOCK2.X series)
- [The HADDOCK2.4 web server: A leap forward in integrative modelling of biomolecular complexes. Nature Prot. 2024](https://www.nature.com/articles/s41596-024-01011-0)

## Development

For development setup and guidelines, see [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md)

## Cite us

If you used `haddock3` for your research, please cite us:

- **Research article**: M. Giulini, V. Reys, J.M.C. Teixeira, B. Jiménez-García, R.V. Honorato, A. Kravchenko, X. Xu, R. Versini, A. Engel, S. Verhoeven, A.M.J.J. Bonvin, [*HADDOCK3: A modular and versatile platform for integrative modelling of biomolecular complexes*](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00969) Journal of Chemical Information and Modeling (2025). doi: 10.1021/acs.jcim.5c00969 [[BioRxiv]](https://www.biorxiv.org/content/10.1101/2025.04.30.651432v1)
 
- **Cite this repository**: M.C. Teixeira, J., Vargas Honorato, R., Giulini, M., Bonvin, A., SarahAlidoost, Reys, V., Jimenez, B., Schulte, D., van Noort, C., Verhoeven, S., Vreede, B., SSchott, & Tsai, R. (2024). haddocking/haddock3: v3.0.0-beta.5 (Version 3.0.0-beta.5) [Computer software]. [https://doi.org/10.5281/zenodo.10527751](https://doi.org/10.5281/zenodo.10527751)
