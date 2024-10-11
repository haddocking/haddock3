"""
Compare CAPRI tables between run folders.

Compare the values in the CAPRI tables looking for missing columns,
missing rows, and different values.

The main function of the script is to compare the integration tests
between two `examples/` folders.

However, this script can also be used to compare two run folders directly.

For more information read 'docs/integration_tests.md'.

## USAGE

To obtain help:

`python compare_runs.py -h`

To compare two `examples/` folders:

`python compare_runs.py -d <PATH-TO-TESTS> -r <PATH-TO-REFERENCE>`

For example, if you are inside the development repository's examples folder,
you need to give `-d` the current path:

`python compare_runs.py -d . -r /home/user/haddock3main/examples`

To compare two specific run folders:

`python compare_runs.py -r1 <PATH-TO-RUN1> -r2 <PATH-TO-RUN2>`

For example:

`python compare_runs.py -r1 scoring/run1-test -r2 scoring/run2-test`

Paramters `-d` and `-r` and mutually exclusive from `-r1` and `r2`.
"""
import argparse
import csv
import os
import sys
from functools import partial
from math import isclose, isnan
from pathlib import Path


try:
    from haddock.gear.config import load as read_config
except Exception:
    print(  # noqa: T201
        "Haddock3 could not be imported. "
        "Please activate the haddock3 python environment.",
        file=sys.stderr,
        )
    sys.exit(1)


printf = partial(print, flush=True)  # noqa: T202


examples = (
    ("docking-antibody-antigen"    , "docking-antibody-antigen-ranairCDR-test.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-ranairCDR-clt-test.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-accessible-test.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-accessible-clt-test.cfg"),  # noqa: E203, E501
    ("docking-antibody-antigen"    , "docking-antibody-antigen-CDR-NMR-CSP-test.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-test.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-mdref-test.cfg"),  # noqa: E203, E501
    ("docking-protein-DNA"         , "docking-protein-DNA-cmrest-test.cfg"),  # noqa: E203, E501
    ("docking-protein-homotrimer"  , "docking-protein-homotrimer-test.cfg"),  # noqa: E203, E501
    ("docking-protein-glycan"      , "docking-protein-glycan-test.cfg"),  # noqa: E203, E501
    ("docking-protein-glycan"      , "docking-protein-glycan-ilrmsd-test.cfg"),  # noqa: E203, E501
    ("docking-protein-glycan"      , "docking-flexref-protein-glycan-test.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand-shape", "docking-protein-ligand-shape-test.cfg"),  # noqa: E203, E501
    ("docking-protein-ligand"      , "docking-protein-ligand-test.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-test.cfg"),  # noqa: E203, E501
    ("docking-protein-peptide"     , "docking-protein-peptide-mdref-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-cltsel-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-protein-protein-mdref-test.cfg"),  # noqa: E203, E501
    ("docking-multiple-ambig"      , "docking-multiple-tbls-test.cfg"),  # noqa: E203, E501
    ("docking-protein-protein"     , "docking-exit-test.cfg"),  # noqa: E203, E501
    ("refine-complex"              , "refine-complex-test.cfg"),  # noqa: E203, E501
    ("peptide-cyclisation"         , "cyclise-peptide-test.cfg"),  # noqa: E203, E501
    ("scoring"                     , "emscoring-test.cfg"),  # noqa: E203, E501
    ("scoring"                     , "mdscoring-test.cfg"),  # noqa: E203, E501
    ("scoring"                     , "emscoring-mdscoring-test.cfg"),  # noqa: E203, E501
    ("analysis"                    , "topoaa-caprieval-test.cfg"),  # noqa: E203, E501
#    ("analysis"                    , "topoaa-clustfcc-test.cfg"),  # noqa: E203, E501
    ("analysis"                    , "topoaa-ilrmsdmatrix-clustrmsd-test.cfg"),  # noqa: E203, E501
    ("analysis"                    , "alascan-test.cfg"),  # noqa: E203, E501
    ("analysis"                    , "contmap-test.cfg"),  # noqa: E203, E501
    )


ap = argparse.ArgumentParser(
    description=__doc__,
    usage="compare_runs.py [-h] [[-d DEVEL] [-r REFERENCE]] OR [[-r1 RUN1] [-r2 RUN2]]",  # noqa: E501
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=(
        "Note: parameters `d` and `-r` are mutually exclusive with "
        "`-r1` and `r2`. You can run either the first case or the "
        "second."
        ),
    )

examples_group = ap.add_argument_group(
    title="Case 1: Compare integration tests in two `examples/` folders",
    )

examples_group.add_argument(
    '-d',
    '--devel',
    type=Path,
    help='Path to the development \'examples/\' folder.',
    default=False,
    )

examples_group.add_argument(
    '-r',
    '--reference',
    type=Path,
    help='Path to the reference \'examples/\' folder.',
    default=False,
    )

specific_run = ap.add_argument_group(
    title="Case 2: Compare two specific run directories",
    )

specific_run.add_argument(
    '-r1',
    '--run1',
    type=Path,
    help='Path to the target run directory.',
    default=False,
    )

specific_run.add_argument(
    '-r2',
    '--run2',
    type=Path,
    help='Path to the reference run directory.',
    default=False,
    )


def load_args():
    """Load argparse arguments."""
    args = ap.parse_args()
    if (args.devel or args.reference) \
            and (args.run1 or args.run2):
        raise ap.error(
            "Case 1 and Case 2 are mutually exclusive. "
            "Give `-d` and `-r` OR `-r1` and `-r2`."
            )
    return args


def print_path_list(pathl):
    """Make printy print for list of paths."""
    _ = (str(Path(*p.parts[-2:])) for p in pathl)
    return "    - " + f"{os.linesep}    - ".join(_)


def print_list(llist):
    """Make nice list."""
    return ", ".join(llist)


def compare_capris(dev_run_dir, ref_run_dir, run_name):
    """Compare CAPRI files in run directory."""
    found_errors = False
    dev_capri_folders = sorted(dev_run_dir.glob("*_caprieval"))
    ref_capri_folders = sorted(ref_run_dir.glob("*_caprieval"))

    if [p.name for p in dev_capri_folders] \
            != [p.name for p in ref_capri_folders]:
        msg = (
            "ERROR FOUND: "
            f"caprieval folders differ for {run_name}:{os.linesep}"
            f"    development:{os.linesep}{print_path_list(dev_capri_folders)}{os.linesep}"  # noqa: E501
            f"    reference:{os.linesep}{print_path_list(ref_capri_folders)}{os.linesep}"  # noqa: E501
            )

        printf(msg)
        return

    for dcp, rcp in zip(dev_capri_folders, ref_capri_folders):
        dev_capri_files = sorted(dcp.glob("*.tsv"))
        ref_capri_files = sorted(rcp.glob("*.tsv"))

        if [p.name for p in dev_capri_files] \
                != [p.name for p in ref_capri_files]:
            msg = (
                "ERROR FOUND: "
                f"caprieval files differ for {run_name} "
                f"folder {dcp.name}:{os.linesep}"
                f"    development:{os.linesep}{print_path_list(dev_capri_files)}{os.linesep}"  # noqa: E501
                f"    reference:{os.linesep}{print_path_list(ref_capri_files)}{os.linesep}"  # noqa: E501
                )

            printf(msg)
            return

        for dcf, rcf in zip(dev_capri_files, ref_capri_files):

            with (open(dcf) as fin1, open(rcf) as fin2):
                devel_tsv = csv.reader(fin1, delimiter="\t")
                ref_tsv = csv.reader(fin2, delimiter="\t")

                devel_table = read_tsv(devel_tsv)
                ref_table = read_tsv(ref_tsv)

                error = compare_tables(devel_table, ref_table)
                if error[0]:
                    printf(f"ERROR FOUND comparing: {str(Path(*dcf.parts[-2:]))}")  # noqa: E501
                    found_errors = True
                print_error[error[0]](*error[1:])

    if not found_errors:
        printf("No errors found - OKAY!")


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
        return round(float(v), 3)
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
    try:
        headers = t1["model"]
    except KeyError:
        headers = t1["cluster_rank"]
    key_values = list(t1.keys())[1:]
    for k in key_values:
        for h, v1, v2 in zip(headers, t1[k], t2[k]):
            if isinstance(v1, str):
                if v1 != v2:
                    return (3, k, h, v1, v2)
            elif isnan(v1) or isnan(v2):
                if not (isnan(v1) and isnan(v2)):
                    return (3, k, h, v1, v2)
            elif isinstance(v1, float):
                if not isclose(v1, v2, abs_tol=0.0011):
                    return (3, k, h, v1, v2)

    return (0, )


def error_1(l1, l2):
    """Print error message for row names."""
    is_in_l1 = set(l1).difference(set(l2))
    is_in_l2 = set(l2).difference(set(l1))

    if not is_in_l1 and not is_in_l2:
        raise AssertionError(
            "BUG FOUND: at least one of these sets should have values."
            )

    printf("Keys in capri files differ:")

    if is_in_l1:

        msg = (
            "These keys are only in the development branch: "
            f"{print_list(is_in_l1)}{os.linesep}"
            )

    if is_in_l2:

        msg = (
            "These keys are only in the reference branch: "
            f"{print_list(is_in_l1)}{os.linesep}"
            )

    printf(msg, flush=True)


def error_2(k, l1, l2):
    """Print error message for different number of values."""
    msg = (
        f"The number of values for {k!r} differ:{os.linesep}"
        f"Development: {len(l1)}{os.linesep}"
        f"Reference: {len(l2)}{os.linesep}"
        )
    printf(msg, flush=True)


def error_3(k, h, v1, v2):
    """Print error message for different values."""
    msg = (
        f"These values differ for {h!r} in model {k!r}:{os.linesep}"
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


def compare_examples(examples, devel, reference):
    """Run all the examples."""
    for folder, file_ in examples:

        params = read_config(Path(devel, folder, file_))
        run_dir = params["final_cfg"]["run_dir"]
        run_name = str(Path(folder, run_dir))

        print(  # noqa: T201
            os.linesep,
            f" {run_name} ".center(80, "*"),
            os.linesep,
            flush=True,
            )  # noqa: T201

        dev_run_dir = Path(devel, folder, run_dir)
        ref_run_dir = Path(reference, folder, run_dir)

        compare_capris(dev_run_dir, ref_run_dir, run_name)

    return


def compare_runs(run1, run2):
    """Compare specific runs."""
    run1 = Path(run1).resolve()
    run2 = Path(run2).resolve()
    run_name = run1.name + " vs. " + run2.name

    printf()
    printf("** Comparing HADDOCK3 runs: ")
    printf(f"Target: {str(run1)}")
    printf(f"Reference: {str(run2)}")
    printf()

    compare_capris(run1, run2, run_name)


def main(devel=False, reference=False, run1=False, run2=False):
    """Execute script logic."""
    if devel and reference:
        compare_examples(examples, devel, reference)
        compare_runs(
            Path(devel, "docking-protein-protein", "run2"),
            Path(reference, "docking-protein-protein", "run2"),
            )
    elif run1 and run2:
        compare_runs(run1, run2)
    else:
        raise ValueError(
            "Parameters are incorrect. "
            "Provide `devel` AND `reference` OR `run1` AND `run2`."
            )


if __name__ == "__main__":
    cmd = load_args()
    main(**vars(cmd))
