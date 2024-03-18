#!/usr/bin/env python3
"""
Copy steps to a new run.

In HADDOCK3 you can copy successful steps from a run directory to a new
directory and use them as starting points for a new run.

Considering the example::

    run1/
        0_topoaa/
        1_rigidbody/
        2_caprieval/
        3_seletop/
        4_flexref/
        (etc...)

You can use `4_flexref` step folder as a starting point for a new run.

USAGE::

    haddock3-copy -r <run_dir> -m <num_modules> -o <new_run_dir>
    haddock3-copy -r run1 -m 0 4 -o run2

Where, ``-m 0 4`` will copy ``0_topoaa`` and ``4_flexref`` to <new_run_dir>.

**Note:** If the new run uses CNS-dependent modules, you **also need**
to copy the folder corresponding to the initial topology creation (the
`topoaa` module).

`haddock3-copy` will also copy the corresponding files in the `data`
directory and update the file contents in the copied folder such that
the information on the run directory and the new step folder names match.
The output result of the above commands is::

    run2/
        0_topoaa/
        1_flexref/

Following, you can use the `haddock3` command with the `--extend-run`
option to continue a new run::

    haddock3 new-config.cfg --extend-run run2
"""
import argparse
import os
import sys

from haddock import log
from haddock.core.typing import ArgumentParser, Callable, FilePath, Namespace
from haddock.libs.libcli import add_version_arg


# Command line interface parser
ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    "-r",
    "--run-dir",
    help="The input run directory.",
    required=True,
    )

ap.add_argument(
    "-m",
    "--modules",
    nargs="+",
    help="The IDs of the steps to copy (separated by spaces).",
    required=True,
    type=int,
    )

ap.add_argument(
    "-o",
    "--output",
    help="The new run directory.",
    required=True,
    )

add_version_arg(ap)


def _ap() -> ArgumentParser:
    return ap


def load_args(ap: ArgumentParser) -> Namespace:
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap: ArgumentParser, main: Callable[..., None]) -> None:
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli() -> None:
    """Execute main client."""
    cli(ap, main)


def main(run_dir: FilePath, modules: list[int], output: FilePath) -> None:
    """
    Copy steps from a run directory to a new run directory.

    Also updates paths references in step files accordingly.

    Parameters
    ----------
    run_dir : str or Path
        Path to the original run directory.

    modules : list of ints
        List of the integer prefix of the modules to copy.

    output : str or Path
        The new run directory to create and where to copy the steps.
    """
    from pathlib import Path

    from haddock.gear.clean_steps import unpack_compressed_and_archived_files
    from haddock.gear.extend_run import (
        copy_renum_step_folders,
        update_contents_of_new_steps,
        )
    from haddock.gear.zerofill import zero_fill
    from haddock.modules import get_module_steps_folders

    log.info("Reading input run directory")
    # get the module folders from the run_dir input
    steps = get_module_steps_folders(run_dir)
    selected_steps = [steps[i] for i in range(len(steps)) if i in modules]
    log.info(f"selected steps: {', '.join(selected_steps)}")

    # make new run dir
    outdir = Path(output)
    try:
        outdir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        log.error(f"Directory {str(outdir.resolve())} already exists.")
        sys.exit(1)
    log.info(f"Created directory: {str(outdir.resolve())}")

    # copy folders over
    zero_fill.set_zerofill_number(len(selected_steps))
    new_step_folder = copy_renum_step_folders(run_dir, outdir, selected_steps)

    # copy data folders
    # `data_steps` are selected to avoid FileNotFoundError because some steps
    # do not have a `data` folder.
    # See https://github.com/haddocking/haddock3/issues/559
    data_steps = [
        step for step in selected_steps
        if os.path.exists(Path(run_dir, "data", step))
        ]
    copy_renum_step_folders(
        Path(run_dir, "data"),
        Path(outdir, "data"),
        data_steps,
        )

    # update step names in files
    # update run dir names in files
    unpack_compressed_and_archived_files(
        new_step_folder,
        ncores=1,
        dec_all=True,
        )
    update_contents_of_new_steps(selected_steps, run_dir, outdir)

    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
