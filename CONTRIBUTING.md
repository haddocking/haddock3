# Contributing to HADDOCK3

Welcome, we made many efforts to facilitate your contribution to this
fantastic project. There are several ways to contribute:

* you can improve tutorials and/or documentation
* improve the code itself (maybe you even found some bug :bug:?)
* improve error messages so they become clearer
* add a new simulation module altogether
* write more unittests (we dare you to do that :godmode:)

Before attempting any development, please install HADDOCK3 following the
instructions in the [INSTALL](INSTALL.md) file.

## Contributing with new code

HADDOCK3 has two main testing workflows. Here, at the repository, we
test the HADDOCK3's Python shell, code style, and building.  Our
Continuous Integration (CI) pipeline is based on
[tox](https://tox.wiki/en/latest/index.html) and GitHub Actions.

To contribute to the HADDOCK3's Python shell, follow these steps:

1.  [Fork](https://www.earthdatascience.org/workshops/intro-version-control-git/about-forks/) the HADDOCK3 repository
2.  Create a new branch in your fork (`git checkout -b <new_branch_name>`)
3.  Develop your code and tests:

    1.  HADDOCK3 source is in `src/haddock`. Always implement code in
    the lowest Python version supported (this case is 3.9)
    2.  Tests sit in the `tests/` folder. Use
    [pytest](https://docs.pytest.org/en/6.2.x/) for testing.

4.  While you are developing (or when you think you have finished), you can
    (should) use our `tox` environments to test your code. Use the
    following commands from the main repository folder:

    1.  `tox -e py39` runs tests in Python 3.9 environment.
    2.  `tox -e lint` shows you errors in the code style.
    3.  `tox -e build` simulates building the HADDOCK3 package.
    4.  If you want to submit high-quality code, you may wish to run
    `tox -e radon` to assess your code cyclomatic complexity.

5.  You can work on these `tox` tests until they all pass green before
submitting your PR.
6.  We also have an `examples` folder with `test` cases that you can run
(should) to ensure the integrity of the python shell as a whole:
    1.  Navigate to any of the examples folder and run the `-test.cfg`
    file (see [USAGE](USAGE.md)).
    2.  If you have a powerful computer and want to run all tests in a
    row, navigate to the `examples` folder and run `python run_tests.py
    -b`. The `-b` flag ensures the run will stop asap an error is found.

7.  Add a list of your new additions to the `CHANGELOG.md` file by
adding a new sub-header as described bellow. This is mandatory for `tox
-e build` to pass.

```markdown
# Changelog

## new_version

* your change one
* your change two

(... the rest of the file ...)
```

7.  If you have difficulties with `tox`, let us know. These `tox` tests
are the same that run on the GitHub Actions once you send the PR. So,
sending the PR is another way to ensure all tests pass.

8.  Submit your PR if you haven't done so yet :wink:

    1.  Now, GitHub allows submitting "Draft PRs" you can use this
    option to let us know you are working on some new good stuff so we
    can help you from the beginning.

## Contributing with documentation

HADDOCK3 has (will have) two sources of documentation: 1) the library
documentation itself that is built from code docstrings and `.rst` pages
in the `docs` folder and 2) the tutorials page at
<https://bonvinlab.org>. If you want to suggest a new tutorial, please
talk with us before start working. On the other hand, if you want to
improve code's documentation, go for it following the same procedure as
described above for *code contributions*.

## You should also know that...

HADDOCK3 folders, files, and code follow a well defined structure with
very specific patterns. Inside each source folder you will (likely) find
a `README` file describing the structure of the folder and presenting
guidelines on how to better contribute to that folder and respective
files.
