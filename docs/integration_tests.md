# Running integration tests manually

*This documentation is for developers.*

The integration tests ensure HADDOCK3 workflows work correctly by performing
short docking or analysis on testing systems. With these tests, we evaluate if
the final scoring values in the developing branch are equal to those obtained in
the `main` branch.

You can run the integration tests locally using our helper scripts and
following the guidelines here.

For a pull request to be accepted, it must pass these integration tests. So we
kindly ask you to run them before submitting a pull request. We are currently
working on incorporating these tests into our automatic CI pipelines.

Running the integration tests for a single branch takes about 1 hour using 8 CPU
cores in a AMD EPYC 7451 24-Core Processor. In case you do not have the computer
power to run these tests, ask us and we will run them for you.

**How to correctly set up and run the integration tests?**

The easiest way is to have two installations of haddock3:
one fixed at the `main` branch and another with your developments.

## 1. Install HADDOCK3 main branch separately

First, you should install HADDOCK3 main branch in a separate folder and
python environment:

```bash
git clone --recursive https://github.com/haddocking/haddock3 haddock3main
cd haddock3main
```

Rename the conda env name in requirements.yml file

```
sed -i 's/haddock3/haddock3main/g' requirements.yml
```

Proceed to install the HADDOCK3 main following the installation
instructions in the repository.

This repository and environment will be the ones you will use to run the
test examples from the `main` branch.

Once installed, you don't need to install this "main" repository again for
future occasions. You can keep using it to test all your developments.

## 3. Running the tests in the `main`

Before running the integration tests, ensure the reference repository `main`
branch is up-to-date.

```bash
cd haddock3main
git checkout main
git pull
```

Remember to activate the corresponding python environment.

```bash
conda activate haddock3main
```

Or, if you used `venv`:

```bash
source venv/bin/activate
```

Navigate to the `examples/` folder and execute the tests:

```bash
cd examples/
python run_tests.py -b
```

This script will run all the test cases in the `examples/` folder (`*-test.cfg`
files) and stop if any of the tests break. For help, see `python run_tests.py
-h`.

## 4. Running the tests on your branch

Navigate to your development repository and activate the haddock3 environment -
you are likely on a different terminal window.

Run the tests the same way as before by navigating to the `examples/` folder and
executing the tests:

```bash
cd examples/
python run_tests.py -b
```

## 5. Comparing both runs

Once tests in the `main` and your development branch complete, that is, all
`*-test.cfg` files were executed, you can compare their CAPRI scores.

Inside the `examples/` folder of your development branch, run:

```bash
python compare-runs.py -r <path-to-haddock-main-examples-folder>
```

This command will tell you if there are any differences in the CAPRI tables:
missing rows or columns and different values. Numeric scores are evaluated with
0.001 tolerance.

If there are differences, likely something went wrong unless you are developing
code that explicitly affects the CAPRI scores. If that is the case, explain it
in the pull request description.

## Additional considerations

### Tests breaking with `run_tests.py`

Tests break during the `run_test.py` step if the new code breaks the core
HADDOCK3 workflow. You should work towards correcting those errors first. Once
the workflows execute completely, you can then compare the CAPRI tables.

### Running on a cluster

Following the instructions here, the test cases will run locally on your
computer. But if you have access to a job scheduling system in a cluster, you
can prepare a `job` script and launch the job as you would to for any other
related process.

### Running individual test cases

If you want to compare a particular test case only, for example, two different
runs for the `docking-protein-protein/docking-protein-protein-test.cfg`, run:

```bash
python compare_runs.py -r1 <PATH-TO-RUN1> -r2 <PATH-TO-RUN2>
```

*Thanks, and enjoy developing and using HADDOCK3 :-)*
