"""
Compare CAPRI tables in two example/ folders.

For more information read 'docs/integration_tests.md'.

USAGE:
    $ python compare_runs.py -h
    $ python compare_runs.py -r <path-to-reference-examples-folder>
"""
import argparse
import csv
import os
import sys
from functools import partial
from math import isclose
from pathlib import Path


try:
    from haddock.gear.config_reader import read_config
except Exception:
    print(  # noqa: T001
        "Haddock3 could not be imported. "
        "Please activate the haddock3 python environment.",
        file=sys.stderr,
        )
    sys.exit(1)


printf = partial(print, flush=True)  # noqa: T002


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
    )


ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-d',
    '--devel',
    type=Path,
    help='Path to the development \'examples/\' folder.',
    default=Path.cwd(),
    )

ap.add_argument(
    '-r',
    '--reference',
    type=Path,
    help='Path to the reference \'examples/\' folder.',
    required=True,
    )


def load_args():
    """Load argparse arguments."""
    return ap.parse_args()


def compare_capris(dev_run_dir, ref_run_dir, run_name):
    """Compare CAPRI files in run directory."""
    dev_capri_folders = sorted(dev_run_dir.glob("*_caprieval"))
    ref_capri_folders = sorted(ref_run_dir.glob("*_caprieval"))

    if dev_capri_folders != ref_capri_folders:
        msg = (
            f"caprieval folders differ for {run_name}:{os.linesep}"
            f"development: {dev_capri_folders}{os.linesep}"
            f"reference: {ref_capri_folders}{os.linesep}"
            )

        printf(msg)
        return

    for dcp, rcp in zip(dev_capri_folders, ref_capri_folders):
        dev_capri_files = sorted(dcp.glob("*.tsv"))
        ref_capri_files = sorted(rcp.glob("*.tsv"))

        if dev_capri_files != ref_capri_files:
            msg = (
                f"caprieval files differ for {run_name} "
                f"folder {dcp.name}:{os.linesep}"
                f"development: {dev_capri_files}{os.linesep}"
                f"reference: {ref_capri_files}{os.linesep}"
                )

            printf(msg)
            return

        for dcf, rcf in zip(dev_capri_files, ref_capri_files):

            with open(dcf) as fin:
                devel_tsv = csv.reader(fin, delimiter="\t")

            with open(rcf) as fin:
                ref_tsv = csv.reader(fin, delimiter="\t")

            devel_table = read_tsv(devel_tsv)
            ref_table = read_tsv(ref_tsv)

            error = compare_tables(devel_table, ref_table)
            print_error[error[0]](error[1:])


def read_tsv(tsv):
    """Read the tsv to a dictionary."""
    d = {}
    for line in tsv:
        values = [try_float(v) for v in line[1:]]
        d[line[0]] = values
    return d


def try_float(v):
    """Try converting string to float."""
    try:
        return float(v)
    except ValueError:
        return v


def compare_tables(t1, t2):
    """Compare two capri tables."""
    # compares if row names are equal
    if list(t1.keys()) != list(t2.keys()):
        return (1, list(t1.keys()), list(t2.keys()))

    # compare if rows (models) have the same number of values
    for k in t1.keys():
        if len(t1[k]) != len(t2[k]):
            return (2, k, t1[k], t2[k])

    # compare value by value, per row
    headers = list(t1["model"].values())
    key_values = list(t1.keys())[1:]
    for k in key_values:
        for h, v1, v2 in zip(headers, t1[k], t2[k]):
            if isinstance(v1, str):
                if v1 != v2:
                    return (3, k, h, v1, v2)
            elif isinstance(v1, float):
                if not isclose(v1, v2, abs_tol=0.001):
                    return (3, k, h, v1, v2)

    return (0, )


def error_1(l1, l2):
    """Print error message for row names."""
    msg = (
        "Keys in both capri files differ:{os.linesep}"
        f"Development: {l1}{os.linesep}"
        f"Reference: {l2}{os.linesep}"
        )
    printf(msg, flush=True)


def error_2(k, l1, l2):
    """Print error message for different number of values."""
    msg = (
        f"The number of values for {k} differ:{os.linesep}"
        f"Development: {len(l1)}{os.linesep}"
        f"Reference: {len(l2)}{os.linesep}"
        )
    printf(msg, flush=True)


def error_3(k, h, v1, v2):
    """Print error message for different values."""
    msg = (
        f"These values differ for {h} in model {k}:{os.linesep}"
        f"Development: {v1}{os.linesep}"
        f"Reference: {v2}{os.linesep}"
        )
    printf(msg, flush=True)


def do_nothing(*args, **kwargs):
    """Do nothing."""
    return


print_error = {
    0: do_nothing,
    1: error_1,
    2: error_2,
    3: error_3,
    }


def main(examples, devel, reference):
    """Run all the examples."""
    for folder, file_ in examples:

        run_name = str(Path(folder, file_))

        print(  # noqa: T001
            os.linesep,
            f" {run_name} ".center(80, "*"),
            os.linesep,
            flush=True,
            )  # noqa: T001

        params = read_config(Path(devel, folder, file_))
        dev_run_dir = Path(devel, folder, params["run_dir"])
        ref_run_dir = Path(reference, folder, params["run_dir"])

        compare_capris(dev_run_dir, ref_run_dir, run_name)

    return


if __name__ == "__main__":
    cmd = load_args()
    main(examples, **vars(cmd))
