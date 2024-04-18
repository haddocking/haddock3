#!/usr/bin/env python3
"""
Clean the output of an HADDOCK3 run directory.

The clean process performs file archiving and file compressing
operations.

All `.inp` and `.out` files are deleted except for the first one, which
is compressed to `.gz`. On the other hand, all `.seed` and `.con` files
are compressed and archived to `.tgz` files. Finally, `.pdb` and `.psf`
files are compressed to `.gz`.

The <run_directory> can either be a whole HADDOCK3 run folder or a
specific folder of the workflow step. <ncores> defines the number of
threads to use; by default uses a single core.

Usage::

    haddock3-clean -h
    haddock3-clean <run_directory>
    haddock3-clean run1
    haddock3-clean run1/1_rigidbody
    haddock3-clean run1 -n  # uses all cores
    haddock3-clean run1 -n 2  # uses 2 cores
"""
import argparse
import sys

from haddock.core.typing import (
    ArgumentParser,
    Callable,
    FilePath,
    Namespace,
    Optional,
)
from haddock.libs import libcli


# Command line interface parser
ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

libcli.add_rundir_arg(ap)
libcli.add_ncores_arg(ap)
libcli.add_version_arg(ap)


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


def main(run_dir: FilePath, ncores: Optional[int] = 1) -> None:
    """
    Clean a HADDOCK3 directory.

    Usually, this concerns compressing and archiving files (see below).

    Parameters
    ----------
    run_dir : str or :external:py:class:`pathlib.Path`.
        The path to the run directory or to a folder of a specific step
        of the workflow.

    ncores : int, or None
        The number of cores to use. If ``None``, use all possible threads.
        Defaults to 1.

    See Also
    --------
    `haddock.gear.clean_steps`
    """
    # anti-pattern to speed up CLI initiation
    from pathlib import Path

    from haddock import log
    from haddock.gear.clean_steps import clean_output
    from haddock.libs.libtimer import log_time
    from haddock.libs.libutil import parse_ncores
    from haddock.modules import get_module_steps_folders, is_step_folder

    log.info(f"Compressing {str(run_dir)!r} folder")
    ncores = parse_ncores(ncores)

    if is_step_folder(run_dir):
        with log_time("compressing took"):
            clean_output(run_dir, ncores)

    else:
        step_folders = get_module_steps_folders(run_dir)
        for folder in step_folders:
            with log_time("compressing took"):
                clean_output(Path(run_dir, folder), ncores)

    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
