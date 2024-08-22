# `haddock3`

<!-- <p align="center">

  <img src="docs/figs/HADDOCK3-logo.png" alt="haddock3-logo" style="vertical-align: middle;">
</p>
<src > -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10527751.svg)](https://doi.org/10.5281/zenodo.10527751)
[![Research Software Directory](https://img.shields.io/badge/rsd-haddock3-00a3e3.svg)](https://research-software-directory.org/software/haddock3)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/8844/badge)](https://www.bestpractices.dev/projects/8844)

[![unit tests](https://github.com/haddocking/haddock3/workflows/tests/badge.svg?branch=main)](https://github.com/haddocking/haddock3/actions?workflow=tests)
[![build](https://github.com/haddocking/haddock3/workflows/build/badge.svg?branch=main)](https://github.com/haddocking/haddock3/actions?workflow=build)
[![docs](https://github.com/haddocking/haddock3/workflows/pages/badge.svg?branch=main)](https://github.com/haddocking/haddock3/actions?workflow=pages)

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e11e7f45400f4e8589cdf5941f95233a)](https://app.codacy.com/gh/haddocking/haddock3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e11e7f45400f4e8589cdf5941f95233a)](https://app.codacy.com/gh/haddocking/haddock3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

![haddock3-logo](docs/figs/HADDOCK3-logo.png)

## Introduction

HADDOCK stands for **H**igh **A**mbiguity **D**riven protein-protein **DOCK**ing, it is a widely used computational tool for modeling protein-protein interactions. Developed by researchers at [Utrecht University](https://uu.nl) in the [BonvinLab](https://bonvinlab.org) it integrates various types of experimental data to guide the docking process, including biochemical, biophysical, and bioinformatics information.

## Installation

Please check the [INSTALL](docs/INSTALL.md) file for instructions.

You might also want to check the following utilities:

- [`haddock-restraints`](https://github.com/haddocking/haddock-restraints): Tool to generate restraints to be used in `haddock3`.
- [`haddock-runner`](https://github.com/haddocking/haddock-runner): Tool to run large scale `haddock3` simulations using multiple input molecules in different scenarios
- [`haddock-tools`](https://github.com/haddocking/haddock-tools): Set of useful utility scripts developed over the years by the BonvinLab group members

## Usage

The most basic usage is;

```bash
haddock3 <configuration-file.toml>
```

Check the [EXAMPLES](examples/README.md) page for more some usage examples and the [MANUAL]() (_coming soon!_) for a more detailed explanation of the configuration file.

## Support

If you encounter any code-related issues, [please open an issue](https://github.com/haddocking/haddock3/issues/new/choose).

If you have any other questions or need help, please contact us at [ask.bioexcel.eu](https://ask.bioexcel.eu/).

If you clone this repository and use `haddock3` for your research, please support us by signing up in [this form](https://forms.gle/LCUHiYHh1hE9rd8L6). This will allow us contact you when needed for `haddock3`-related issues, and also provide us a mean to demonstrate impact when reporting for grants - which grealty helps us to keep the project alive!

## Development

🚧 _Coming soon! We are currently working on the development guide_ 🚧

### Code Documentation

The code documentation is automatically built and hosted at [bonvinlab.org/haddock3](https://www.bonvinlab.org/haddock3/).

To build it locally (considering you have followed the installation instructions):

```bash
tox -e docs
```

> Warning messages are expected, but the documentation should be built successfully.

The rendered documentation will be available at `haddock3-docs/index.html`. This will open a local webpage with the

### Contributing

Check the [CONTRIBUTING](CONTRIBUTING.md) file for instructions on how to contribute with the project!

<!-- ---

Happy HADDOCking!

<img src="https://www.bonvinlab.org/images/bio-haddock.png" alt="haddock" width="50px"> -->
