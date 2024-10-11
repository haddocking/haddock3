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

    $ python run_tests.py -h
    $ python run_tests.py     # runs all examples regardless of errors
    $ python run_tests.py -b  # stops asap an error is found
"""

import argparse
import os
import subprocess
import sys
from functools import partial
from shutil import rmtree


try:
    from haddock.gear.config import load as read_config
    from haddock.libs.libio import working_directory
except Exception:
    print(  # noqa: T201
        "Haddock3 could not be imported. "
        "Please activate the haddock3 python environment.",
        file=sys.stderr,
    )
    sys.exit(1)


# edit this tuple to add or remove examples.
# keys are the examples folder, and values are the configuration files
# the whitespaces below are anti-pythonic but facilitate reading :-)
examples = (
    (
        "docking-antibody-antigen",
        "docking-antibody-antigen-ranairCDR-test.cfg",
    ),  # noqa: E203, E501
    (
        "docking-antibody-antigen",
        "docking-antibody-antigen-ranairCDR-clt-test.cfg",
    ),  # noqa: E203, E501
    (
        "docking-antibody-antigen",
        "docking-antibody-antigen-CDR-accessible-test.cfg",
    ),  # noqa: E203, E501
    (
        "docking-antibody-antigen",
        "docking-antibody-antigen-CDR-accessible-clt-test.cfg",
    ),  # noqa: E203, E501
    (
        "docking-antibody-antigen",
        "docking-antibody-antigen-CDR-NMR-CSP-test.cfg",
    ),  # noqa: E203, E501
    ("docking-protein-DNA", "docking-protein-DNA-test.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA", "docking-protein-DNA-mdref-test.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA", "docking-protein-DNA-cmrest-test.cfg"),
    (
        "docking-protein-homotrimer",
        "docking-protein-homotrimer-test.cfg",
    ),  # noqa: E203, E501
    ("docking-protein-glycan", "docking-protein-glycan-test.cfg"),  # noqa: E203, E501
    (
        "docking-protein-glycan",
        "docking-protein-glycan-ilrmsd-test.cfg",
    ),  # noqa: E203, E501
    (
        "docking-protein-glycan",
        "docking-flexref-protein-glycan-test.cfg",
    ),  # noqa: E203, E501
    (
        "docking-protein-ligand-shape",
        "docking-protein-ligand-shape-test.cfg",
    ),  # noqa: E203, E501
    ("docking-protein-ligand", "docking-protein-ligand-test.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide", "docking-protein-peptide-test.cfg"),  # noqa: E203, E501
    (
        "docking-protein-peptide",
        "docking-protein-peptide-mdref-test.cfg",
    ),  # noqa: E203, E501
    ("docking-protein-protein", "docking-protein-protein-test.cfg"),  # noqa: E203, E501
    (
        "docking-protein-protein",
        "docking-protein-protein-cltsel-test.cfg",
    ),  # noqa: E203, E501
    (
        "docking-protein-protein",
        "docking-protein-protein-mdref-test.cfg",
    ),  # noqa: E203, E501
    ("docking-multiple-ambig", "docking-multiple-tbls-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein", "docking-exit-test.cfg"),  # noqa: E203, E501
    ("refine-complex", "refine-complex-test.cfg"),  # noqa: E203, E501
    ("peptide-cyclisation" , "cyclise-peptide-test.cfg"),  # noqa: E203, E501
    ("scoring", "emscoring-test.cfg"),  # noqa: E203, E501
    ("scoring", "mdscoring-test.cfg"),  # noqa: E203, E501
    ("scoring", "emscoring-mdscoring-test.cfg"),  # noqa: E203, E501
    ("analysis", "topoaa-caprieval-test.cfg"),  # noqa: E203, E501
    ("analysis", "topoaa-ilrmsdmatrix-clustrmsd-test.cfg"),  # noqa: E203, E501
    ("analysis", "alascan-test.cfg"),  # noqa: E203, E501
    ("analysis", "contmap-test.cfg"),  # noqa: E203, E501
)


ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

ap.add_argument(
    "-b",
    "--break-on-errors",
    action="store_true",
    help=(
        "Stop execution as soon an example gives an error. "
        "If not given, runs all examples regardless of errors."
    ),
)


def load_args():
    """Load argparse arguments."""
    return ap.parse_args()


def _run_subprocess_cmd(cmd_: str, check: bool = True) -> None:
    """Run provided command line with subprocess.

    Parameters
    ----------
    cmd_ : str
        The Haddock3 command line
    check : bool, optional
        Should the exit status be checked, by default True
    """
    subprocess.run(
        cmd_,
        shell=True,
        check=check,
        stdout=sys.stdout,
        stderr=sys.stderr,
    )


def main(examples, break_on_errors=True):
    """Run all the examples."""
    # Preset subprocess arguments
    run_subprocess_cmd = partial(_run_subprocess_cmd, check=break_on_errors)

    # Loop over config files
    for folder, file_ in examples:

        print(  # noqa: T201
            os.linesep,
            f" {file_.upper()} ".center(80, "*"),
            os.linesep,
            flush=True,
        )  # noqa: T201

        with working_directory(folder):

            # obtain run directory
            all_params = read_config(file_)
            rundir = all_params["final_cfg"]["run_dir"]
            # remove eventual previous run
            rmtree(rundir, ignore_errors=True)

            # run example
            run_subprocess_cmd(f"haddock3 {file_}")

            # test sub-commands / parameters on the prot-prot run
            if file_ == "docking-protein-protein-test.cfg":
                # perform a restart step
                run_subprocess_cmd(f"haddock3 {file_} --restart 5")

                # perform a restart step from 0
                run_subprocess_cmd(f"haddock3 {file_} --restart 0")

                # test copy run
                rmtree("run2", ignore_errors=True)
                run_subprocess_cmd("haddock3-copy -r run1-test -m 0 4 -o run2")

                # test exit with extend-run
                rmtree("run2", ignore_errors=True)
                run_subprocess_cmd("haddock3-copy -r run1-test -m 0 4 -o run2")
                run_subprocess_cmd(
                    "haddock3 docking-extend-run-exit-test.cfg --extend-run run2",  # noqa: E501
                )

                # test exit with --restart
                rmtree("run1-restart-exit-test", ignore_errors=True)
                run_subprocess_cmd("cp -r run1-test run1-restart-exit-test")
                run_subprocess_cmd(
                    "haddock3 docking-restart-exit-test.cfg --restart 3",
                )

                # Copy run for haddock3-re commands
                rmtree("run1-re", ignore_errors=True)
                run_subprocess_cmd(
                    "haddock3-copy -r run1-test -m 0 7 9 -o run1-re",
                )

                # perform a haddock3 re-scoring command
                run_subprocess_cmd(
                    "haddock3-re score -e 1.1 -w 1 -d 0.3 -b 1 -a 1 run1-re/2_caprieval",  # noqa : E501
                )

                # perform a haddock3 re-clustfcc command
                run_subprocess_cmd(
                    "haddock3-re clustfcc -f 0.5 -s 0.7 -t 2 run1-re/1_clustfcc",  # noqa : E501
                )

                # FIXME: Make this runs properly function
                # perform a haddock3 re-clustrmsd command

                # perform haddock3 --extend-run on re-run
                # run_subprocess_cmd(
                #     "haddock3 docking-re-extend-run-test.cfg --extend-run run1-re",  # noqa : E501
                #     )

                # perform haddock3 --restart on re-run
                # run_subprocess_cmd(
                #     "haddock3 docking-re-restart-test.cfg --restart 5",
                #     )

    return


if __name__ == "__main__":
    cmd = load_args()
    main(examples, **vars(cmd))
