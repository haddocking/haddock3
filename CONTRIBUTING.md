Here at HADDOCK3, we have made it such that it is easy to contribute to
the project. There are several ways to contribute: you can improve
tutorials and documentation, the code itself, add new simulation
modules, etc.

Before attempting any development, please install HADDOCK3 following the
instructions in the INSTALL file.

**Contributing with new code**

HADDOCK3 has two main testing workflows. Here, at the repository, we
test the HADDOCK3's Python shell, code quality and style and building.
Our Continuous Integration (CI) pipeline is based on [tox](https://tox.wiki/en/latest/index.html) and GitHub
Actions. To benchmark modules (using CNS), we have a private Jenkins
server.

To contribute to the HADDOCK3's Python shell, follow these steps:

1.  [Fork](https://www.earthdatascience.org/workshops/intro-version-control-git/about-forks/) the HADDOCK3 repository

2.  Create a new branch in your fork (`git checkout -b <new_branch_name>`)

3.  Develop your code and tests

    1.  HADDOCK3 source is in `src/haddock`. Always implement code in the lowest Python version supported (this case is 3.8)

    2.  Tests sit in the `tests/` folder. Use [pytest](https://docs.pytest.org/en/6.2.x/) for testing.

4.  While you are developing (or when you think you finished), you can
    (should) use our `tox` environments to test your code. Use the
    following commands from the main repository folder:

    1.  `tox -e py38` runs tests in Python 3.8 environment.

    2.  `tox -e lint` shows you errors in the code style.

    3.  `tox -e build` simulates building the HADDOCK3 package.

    4.  `tox -e prreqs` adds further confirmations you don't forget anything before sending your Pull Request.

    5.  If you want to submit high-quality code, you may wish to run `tox -e radon` to assess your code cyclomatic complexity.

5.  You can work on these `tox` tests until they all pass green before submitting your PR.

6.  If you have difficulties with `tox`, let us know. These `tox` tests are the same that run on the GitHub Actions once you send the PR. So, sending the PR is another way to ensure all tests pass.

7.  Submit your PR if you haven't done so yet :wink:

    1.  Now, GitHub allows submitting "Draft PRs" you can use this option to let us know you are working on some new good stuff so we can help you from the beginning.

**Contributing with documentation:**

We currently have not defined what will be the documentation platform
for HADDOCK3. So, stay tuned for that :wink:; meanwhile, you can always
contribute with code documentation if you find it necessary.
