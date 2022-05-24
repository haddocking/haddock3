Welcome HADDPCK3-Beta version

The `main` branch represents the latest state of HADDOCK v3. Currently,
stable beta version.

* * *

# HADDOCK3

## 1. Installation

To install HADDOCK3 follow the instructions in the [INSTALL](docs/INSTALL.md) file.

## 2. Documentation

HADDOCK3-beta documentation is not yet hosted online. You need to generate
it locally. First, install HADDOCK3 and activate the `haddock3` python
environment. Then, in your terminal window, run:

```bash
tox -e docs
```

After, use your favorite browser to open the file
`dist/docs/index.html`. Navigate around, enjoy, and contribute.

## 2. Examples

### 2.1. Basic scoring of an ensemble of 5 structures:

In the `examples/` folder you find several examples for you to test and
learn HADDOCK3. Additional information is in the documentation pages.

```bash
cd examples/recipes/scoring/
haddock3 scoring.toml
```

## 3. Contribute

If you want to contribute to HADDOCK3's development, read the
[CONTRIBUTING](CONTRIBUTING.md) file for instructions.

## 4. Keep in contact and support us

HADDOCK3 is an academic project supported by various grants, including the EU BioExcel Center of Excellence for Computational Biomolecular Research. HADDOCK3 is fully open-source and free to download. If you clone this repository and use HADDOCK3 for your research, please support us by signing [this Goole form][googleform] if you have not yet done so. This will allow us contact you when needed for HADDOCK3-related issues, and also provide us a mean to demonstrate impact when reporting for grants.

[googleform]: https://docs.google.com/forms/d/e/1FAIpQLScDcd0rWtuzJ_4nftkDAHoLVwr1IAVwNJGhbaZdTYZ4vWu25w/viewform
