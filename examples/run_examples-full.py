"""
Run all examples in a row.

This script should be executed from inside the haddock3/examples/ folder.

Inside this script, we have listed only the examples that are supposed to run
locally, meaning examples referring to HPC, or other special cases,  won't be
executed. In other words, examples are not fetch automatically from the folders,
and new examples must be added (or deleted) manually from inside the script.

Run folders that overlap with the `run_dir` parameter in the configuration files
will be deleted.

If you see errors related to python import statements, make sure you have
the haddock3 environment activated.

USAGE:

    $ python run_examples.py -h
    $ python run_examples.py     # runs all examples regardless of errors
    $ python run_examples.py -b  # stops asap an error is found
"""
import argparse
import subprocess
import sys
from shutil import rmtree


try:
    from haddock.libs.libio import working_directory
    from haddock.gear.config_reader import read_config
except Exception:
    print(  # noqa: T001
        "Haddock3 could not be imported. "
        "Please activate the haddock3 python environment.",
        file=sys.stderr,
        )
    sys.exit(1)


# edit this dictionary to add or remove examples.
# keys are the examples folder, and values are the configuration files
# spacings are anti-pythonic but facilitate reading :-)
examples = (
    ("docking-protein-DNA"         , "docking-protein-DNA-full.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-mdref-full.cfg"),  # noqa: E203, E501
    ("docking-protein-homotrimer"  , "docking-protein-homotrimer-full.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand-shape", "docking-protein-ligand-shape-full.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand"      , "docking-protein-ligand-full.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-full.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-mdref-full.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-full.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-mdref-full.cfg"),  # noqa: E203, E501
    )


ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-b',
    '--break-on-errors',
    action='store_true',
    help=(
        "Stop execution as soon an example gives an error. "
        "If not given, runs all examples regardless of errors."
        ),
    )


def load_args():
    """Load argparse arguments."""
    return ap.parse_args()


def main(examples, break_on_errors=True):
    """Run all the examples."""
    for folder, file_ in examples:

        print()  # noqa: T001
        print(f" {file_.upper()} ".center(80, "*"))  # noqa: T001
        print()  # noqa: T001

        with working_directory(folder):

            params = read_config(file_)
            rmtree(params["run_dir"], ignore_errors=True)
            subprocess.run(
                f"haddock3 {file_}",
                shell=True,
                check=break_on_errors,
                stdout=sys.stdout,
                stderr=sys.stderr,
                )

    return


if __name__ == "__main__":
    cmd = load_args()
    main(examples, **vars(cmd))
