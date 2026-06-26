#!/usr/bin/env python3
"""
Unpack the output of an HADDOCK3 run directory.

The unpack process performs file unpacking and file decompressing
operations.  File with extension `seed` and `con` are unpacked from
their `.tgz` files.  While files with `.pdb.gz` and `.psf.gz` extension
are uncompressed.  If `--all` is given, unpack also `.inp.gz` and
`.out.gz` files.

This CLI performs the opposite operations as the ``haddock3-clean``
command-line.

The <run_directory> can either be a whole HADDOCK3 run folder or a
specific folder of the workflow step. <ncores> defined the number of
threads to use.

Usage::

    haddock3-unpack -h
    haddock3-unpack -r <run_directory>
    haddock3-unpack run1
    haddock3-unpack run1/1_rigidbody
    haddock3-unpack run1 -n  # uses all cores
    haddock3-unpack run1 -n 2  # uses 2 cores
    haddock3-unpack run1 -n 2 -a
    haddock3-unpack run1 -n 2 --all
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

ap.add_argument(
    "--all",
    "-a",
    dest="dec_all",
    help="Unpack all files (includes `.inp` and `.out`).",
    action="store_true",
)

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


def main(run_dir: FilePath, ncores: Optional[int] = 1, dec_all: bool = False) -> None:
    """
    Unpack a HADDOCK3 run directory step folders.

    Usually, this concerns uncompressing and unpacking files.

    It performs the oposity function as
    :py:func:`haddock.clis.cli_clean.main`.

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
    from haddock.gear.clean_steps import unpack_compressed_and_archived_files
    from haddock.libs.libtimer import log_time
    from haddock.libs.libutil import parse_ncores
    from haddock.modules import get_module_steps_folders, is_step_folder

    log.info(f"Unpacking {str(run_dir)!r} folder")
    ncores = parse_ncores(ncores)

    if is_step_folder(run_dir):
        with log_time("unpacking took"):
            unpack_compressed_and_archived_files(
                [run_dir],
                ncores,
                dec_all=dec_all,
            )

    else:
        step_folders = [Path(run_dir, p) for p in get_module_steps_folders(run_dir)]
        with log_time("unpacking took"):
            unpack_compressed_and_archived_files(
                step_folders,
                ncores,
                dec_all=dec_all,
            )

    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
