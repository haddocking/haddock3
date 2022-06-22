"""
Run all test examples in a row.

This script should be executed from inside the `examples/` folder.

This script runs the `*-test.cfg` files. These test cases are not fetched
automatically. If you want to add/remove cases you need to edit this script.
The HADDOCK3 team has defined here only those test cases that are part of the
integration tests.

This script will delete existing run directories which names overlap with those
defined in the test configuration files.

If you see errors related to python import statements, make sure you have
the HADDOCK3 environment activated.

A breaking example means something is wrong in the HADDOCK3 core workflow.
You should work towards solving that problem or contact the HADDOCK3 team.

USAGE:

    $ python run_examples.py -h
    $ python run_examples.py     # runs all examples regardless of errors
    $ python run_examples.py -b  # stops asap an error is found
"""
import argparse
import os
import subprocess
import sys
from shutil import rmtree


try:
    from haddock.gear.config_reader import read_config
    from haddock.libs.libio import working_directory
except Exception:
    print(  # noqa: T201
        "Haddock3 could not be imported. "
        "Please activate the haddock3 python environment.",
        file=sys.stderr,
        )
    sys.exit(1)


# edit this dictionary to add or remove examples.
# keys are the examples folder, and values are the configuration files
# the whitespaces below are anti-pythonic but facilitate reading :-)
examples = (
    ("docking-antibody-antigen"    , "docking-antibody-antigen-ranairCDR-test.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-ranairCDR-clt-test.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-accessible-test.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-accessible-clt-test.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-test.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-mdref-test.cfg"),  # noqa: E203, E501
    ("docking-protein-homotrimer"  , "docking-protein-homotrimer-test.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand-shape", "docking-protein-ligand-shape-test.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand"      , "docking-protein-ligand-test.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-test.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-mdref-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-cltsel-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-mdref-test.cfg"),  # noqa: E203, E501
    ("refine-complex"              , "refine-complex-test.cfg"),  # noqa: E203, E501
    ("scoring"                     , "emscoring-test.cfg"),  # noqa: E203, E501
    ("scoring"                     , "mdscoring-test.cfg"),  # noqa: E203, E501
    ("scoring"                     , "emscoring-mdscoring-test.cfg"),  # noqa: E203, E501
    ("analysis"                    , "topoaa-caprieval-test.cfg"),  # noqa: E203, E501
    ("analysis"                    , "topoaa-clustfcc-test.cfg"),  # noqa: E203, E501
    ("analysis"                    , "topoaa-rmsdmatrix-clustrmsd-test.cfg")  # noqa: E203, E501
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

        print(  # noqa: T201
            os.linesep,
            f" {file_.upper()} ".center(80, "*"),
            os.linesep,
            flush=True,
            )  # noqa: T201

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

            # perform a restart step
            if file_ == "docking-protein-protein-test.cfg":
                subprocess.run(
                    f"haddock3 {file_} --restart 5",
                    shell=True,
                    check=break_on_errors,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    )

            # perform a restart step from 0
                subprocess.run(
                    f"haddock3 {file_} --restart 0",
                    shell=True,
                    check=break_on_errors,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    )

                rmtree("run2", ignore_errors=True)
                subprocess.run(
                    "haddock3-copy -r run1-test -m 0 4 -o run2",
                    shell=True,
                    check=break_on_errors,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    )

                subprocess.run(
                    "haddock3 docking-protein-protein-test-start-from-cp.cfg --extend-run run2",  # noqa: E501
                    shell=True,
                    check=break_on_errors,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    )

    return


if __name__ == "__main__":
    cmd = load_args()
    main(examples, **vars(cmd))
