# Running integration tests manually

Meanwhile we implemented automate integration tests, you can run them
manually using our helper scripts.

The integration tests ensure HADDOCK3 workflows are working properly.
You should run the integration tests every time you propose a new pull
request.

## 1. How to properly setup and run the integration tests?

The easiest way is to have two installations of haddock3:
one fixed at the `main` branch and other with your developments.

Ensure the "HADDOCK3 main" repository and your HADDOCK3
development repository are placed in the same parent folder. For
example:

```
some-folder/
    - haddock3/
    - haddock3main/
```

## 2. Install HADDOCK3 main branch in a separate folder and python environment

```bash
git clone --recursive https://github.com/haddocking/haddock3 haddock3main
cd haddock3main

# rename the conda env name in requirements.yml file
sed -i 's/haddock3/haddock3main/g' requirements.yml
```

Proceed installing the HADDOCK3 main as usual.
This repository and environment will be the one you will use to run the
test examples from the main branch.

Once installed you don't need to install this "main" repository again.

## 3. Running the tests in the "main"

Before running the integration tests, ensure the "main" branch
repository is up-to-date.

```bash
git pull
```

Navigate to the `examples/` folder and execute the tests:

```bash
python run_tests.py -b
```

This script will run all the test cases in the `examples/` directory and
stop imminently if any of the tests break. For help, see `python
run_tests.py -h`.

If you are on a Slurm supported system, you can use the `test.job` file.
You might need to edit it to match the paths and environment names.

## 4. Running the tests on your branch

Navigate to your development repository and haddock3 environment - you
are likely on a different terminal window.

Run the tests the same way:

Navigate to the `examples/` folder and execute the tests:

```bash
python run_tests.py -b
```

If you are on a Slurm supported system, you can use the `test.job` file.

## 5. Comparing both runs

Once both runs completed, that is, all `*-test.cfg` files were executed,
you can compare the different results for the CAPRI scoring modules.

Inside the `examples/` folder for your development branch, run:

```bash
./diff_capri.sh
```

This will tell you if there are any differences in the CAPRI scores. If
there are, likely something went wrong, unless you are developing code
that specifically affects the CAPRI scores.

**Note:** You can configure the relative paths between the two haddock
installations in case you installed `haddock main` in a different
path from the specified here.
