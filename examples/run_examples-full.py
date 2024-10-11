"""
Run all full examples in a row.

This script should be executed from inside the haddock3/examples/ folder.

This script runs the `*-full.cfg` files. These test cases are not fetched
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

    $ python run_examples-full.py -h
    $ python run_examples-full.py     # runs all examples regardless of errors
    $ python run_examples-full.py -b  # stops asap an error is found
"""
import argparse
import subprocess
import sys
from shutil import rmtree


try:
    from haddock.gear.config import load as read_config
    from haddock.libs.libio import working_directory
except Exception:
    print(  # noqa: T001
        "Haddock3 could not be imported. "
        "Please activate the haddock3 python environment.",
        file=sys.stderr,
        )
    sys.exit(1)


# edit this tuple to add or remove examples.
# keys are the examples folder, and values are the configuration files
# spacings are anti-pythonic but facilitate reading :-)
examples = (
    ("docking-protein-DNA"         , "docking-protein-DNA-full.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-cltsel-full.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-mdref-full.cfg"),  # noqa: E203, E501
    ("docking-protein-homotrimer"  , "docking-protein-homotrimer-full.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand"      , "docking-protein-ligand-full.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand-shape", "docking-protein-ligand-shape-full.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-full.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-cltsel-full.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-mdref-full.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-full.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-cltsel-full.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-mdref-full.cfg"),  # noqa: E203, E501
    ("docking-multiple-ambig"      , "docking-multiple-tbls-clt-full.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-NMR-CSP-full.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-accessible-full.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-accessible-clt-full.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-ranairCDR-full.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-ranairCDR-clt-full.cfg"),  # noqa: E203, E501
    ("peptide-cyclisation"         , "cyclise-peptide-full.cfg"),  # noqa: E203, E501
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
            # obtain run directory
            all_params = read_config(file_)
            rundir = all_params['final_cfg']["run_dir"]
            # remove eventual previous run
            rmtree(rundir, ignore_errors=True)

            # run haddock with config file
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
