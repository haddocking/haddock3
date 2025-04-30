# Contributing to HADDOCK3

Welcome, we made many efforts to facilitate your contribution to this
fantastic project. There are several ways to contribute:

* you can improve tutorials and/or documentation
* improve the code itself (maybe you even found some bug :bug:?)
* improve error messages so they become clearer
* add a new simulation module altogether
* write more unittests (we dare you to do that :godmode:)

Before attempting any development, please install HADDOCK3 following the
instructions in the [INSTALL](docs/INSTALL.md) file. Afterwards, follow the
instructions in this file.

## 1. Contributing new code

HADDOCK3 has two main testing workflows. Here, within the repository, we
test the HADDOCK3's Python shell, code style, and package building.  Our
Continuous Integration (CI) pipeline is based on
pytest and GitHub Actions. We will explain you how to use them.

To contribute to the HADDOCK3's Python shell, follow these steps:

1.  [Fork][fork] the HADDOCK3 repository
2.  Create a new branch in your fork (`git checkout -b <new_branch_name>`)
3.  Develop your code:

    1.  HADDOCK3 source is in `src/haddock`. Always implement code in
    the lowest Python version supported (this case is 3.9).
    2. See more details on how to contribute with code and tests in the
    subheadings below.

4.  Test your code:
    1.  Unit tests sit in the `tests/` folder. Use
    [pytest](https://docs.pytest.org/en/6.2.x/) for testing by running `pytest tests/`.
    2.  Integration tests are located in the `integration_tests/` folder.
        You can run them with `pytest integration_tests/`.

5.  You can work on these tests until they all pass green before
submitting your PR.

6. Check if your contribution fulfills the requeriments proposed in the
PR template, these are based on the experience of some developers to ensure the long-term
survival of the codebase. Suggestions are always welcomed.

7.  Add a list of your new additions to the `CHANGELOG.md` file by
adding a new sub-header as described bellow.
Note this applies only after we have released the
stable `3.0.0` version.

```markdown
# Changelog

## new_version

* your change one
* your change two

(... the rest of the file ...)
```

8.  If you have difficulties with the `tests`, let us know. These tests
are the same that run on the GitHub Actions once you send the PR. So,
sending the PR is another way to ensure all tests pass.

9.  Submit your PR if you haven't done so yet :wink:

    1.  GitHub allows submitting "Draft PRs". You can use this
    option to let us know you are working on some new good stuff so we
    can help you from the beginning.

### 1.1 Contributing with code (additional details)

**You should also know that...** HADDOCK3 folders, files, and code
follow a well defined structure with very specific patterns. Inside each
source folder you will (likely) find a `README` file describing the
structure of the folder and presenting guidelines on how to better
contribute to that folder and respective files. Summarizing here:

* New command-line clients go in the `clis/` folder. See how the current
  clients here created and use that as a template. Remember to update also the
  `pyproject.toml` file.
* Add any new functions that you foresee are general and
  could serve different places in the code to the `libs/` folder. Find a
  `lib*.py` files that would serve your purpuses, otherwise create a new
  one.
* Any plug-in like functionality, for example *"check if input is
  correct"*, should go into its own python module inside the `gear/`
  folder. Any variables, functions, and classes related only to that
  implementation should go all inside the same module (Python file).
* Any general hard definitions, or physical constants, should go in the
  `core/` folder.
* If you want to implement a new HADDOCK3 module, navigate to the
  `src/haddock/modules/_template_cat` folder and follow the instructions in the `README.md`
  file there. You will see that all folders and files follow a pattern.
* Talk with us before developing any CNS related part.

### 1.2 Contributing with tests (additional details)

Inside the `tests/` folder you will find several `test_*.py` files.
Normally, each file has the tests for each `*.py` file in the source. If
you create new `*.py` files you should create a new test file of the
same name, `test_new_file_name.py`. Aim at 100% test coverage for the
code you have created. Write tests according to [pytest][pytest]. You
can see examples in our `test_*.py` files. You can run the tests using
the `pytest tests/` command explained above. Or, if you want to run the
tests for a singly file use `pytest tests/test_myfile.py`.

### 1.3 Dependencies

HADDOCK3 is highly interconnected with other projects. Many of them use
HADDOCK3 core functionalities. Therefore, we aim to keep HADDOCK3 with
the lowest dependency footprint possible. Avoid adding dependencies when
developing new functionalities. How?

1. Use the Python standard library as much as possible.
1. Numpy is always allowed.
1. If you need to implement heavy calculations, it is best you use
[Numba][numba] instead of C libraries. Talk with us before.
1. You need a small function from a large library. Try to reimplement it
yourself with the Python standard library.
    1. If you can't, talk first with
us by opening an issue.
1. You need a big and complex function or maybe a whole python file from
another project. If licenses are compatible, copy their code to the
HADDOCK3 repository writing all necessary headers to grant authorship
and comply with license requirements.
    1. If licenses aren't compatible, talk to us. We may have an
    alternative.
1. Your new module largely depends on a library and reimplementing or
copying is not an option. Then, consider if we can use that dependency
as a **runtime (thirdparty) dependency** (like `gdock` or `lightdock`) instead of an
installation dependency.
1. Nothing of the above is possible. You **really** need an *install
dependency*. Talk with us before.

### 1.4 Code structure

Write code in the form of small functions, because small functions are
easier to test. Combine small functions to compose larger functions. If
you need to use a global variable in a function, pass it as a default
value of a parameter. Avoid using complex classes, or avoid using
classes at all, unless you need to maintain a state or you really know
what you are doing. Flat is better than nested (though it's harder to
write, but it's easier to maintain and read). Use long variable names if
needed. Write comments that explain why you do stuff, and not how you do
stuff. Use the `TODO:` flag in your comments to note something for the
future. If needed, raise an issue.

### 1.5 Creating a new module

To develop a new HADDOCK3 module follow our guidelines and templates
under `src/haddock/modules/_template_cat/_template_mod/`.

## 2. Contributing with documentation

You may contribute to HADDOCK3 documentation by improving parts where
documentation is lacking or writing the documentation for the new code
you propose. HADDOCK3 documentation is hosted online at
https://bonvinlab.org/haddock3.

HADDOCK3 documentation is rendered with
[Sphinx](https://www.sphinx-doc.org/en/master/) combining markdown
files, restructured text files, and extracting the docstrings in the
source code.

The `docs/` folder contains all the files used by Sphinx to compile the
documentation to HTML files. To incorporate new documentation pages or
update the existing ones, navigate around the `docs/` folder to learn
how we have structured it and add/edit the files you find relevant. You
will see that the structure of folders and files follows the design of
the documentation website.

You can render the documentation locally to inspect the end result
before creating a pull request. To compile the documentation locally:
- install HADDOCK in development mode with the `docs` extra requirements
  (`pip install -e '.[docs]'`)
- activate the `haddock3` python environment inside the haddock3 github
folder
- run the following commands:
```bash
sphinx-apidoc -f -e -o docs/ src/haddock -d 1
sphinx-build -b html docs haddock3-docs
```

We invite you to read through Sphinx-doc
webpage if you want to exploit any advanced feature of Sphinx, but we
already provide examples for virtually any use you may need.

If you need to install any additional library, talk to us first. The
documentation requirements are in the `devtools` folder.

Finally, if you find the need to generate new pages during the HTML
compilation part, you may follow the `devtools/build_defaults_rst.py` as
an example. See also the `docs/conf.py` file `setup(app)` line.

**Troubleshooting:**

1. If you add any new dependency (import statement) in the code, you
  need to add that library name to the `mock` list in the `docs/conf.py`
  file.

[fork]: https://docs.github.com/en/get-started/quickstart/fork-a-repo
[pytest]: https://docs.pytest.org/ "pytest"
[flake]: https://flake8.pycqa.org/en/latest/ "flake8"
[numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[numba]: https://numba.pydata.org/ "Numba"
