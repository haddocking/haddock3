"""
MPI-wrapper for executing pickled CNSJobs created by libmpi.

This was developed for use of lbimpi but it might be useful in some specific
 scenario as a cli.

For more information please refer to the README.md in the examples folder.

Usage::

    haddock3-mpitask -h
    haddock3-mpitask tasks.pkl
"""

import argparse
import pickle
import sys

from haddock.core.typing import (
    AnyT,
    ArgumentParser,
    Callable,
    FilePath,
    Namespace,
)
from haddock.libs.libsubprocess import CNSJob


# Note! ########################################################################################
# It's important to not define `MPI` and `COMM` directly but to use the `get_mpi`
#  function to ensure that the MPI module is only loaded when needed.
# This is important for testing purposes, and this CLI should further be refactored to not use
#  global variables.
MPI = None
COMM = None
# Note! ########################################################################################


def get_mpi():
    """Lazy load MPI and COMM."""
    global MPI, COMM
    if MPI is None:
        try:
            from mpi4py import MPI as _MPI

            MPI = _MPI
            COMM = MPI.COMM_WORLD
        except ImportError as e:
            _msg = (
                f"{e} - To run this cli you must have mpi4py and "
                "OpenMPI installed in the system"
            )
            sys.exit(_msg)
    return MPI, COMM


def split_tasks(task_l: list[AnyT], n: int) -> list[list[AnyT]]:
    """Split tasks into equal lenght chunks."""
    return [task_l[_i::n] for _i in range(n)]


# ========================================================================#
# helper functions to enhance flexibility and modularity of the CLIs

ap = argparse.ArgumentParser(
    prog="MPI-wrapper for executing pickled CNSJobs created by libmpi",
    description=__doc__,
)

ap.add_argument(
    "pickled_tasks",
    help="The input pickled tasks path",
)


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


# ========================================================================#


def main(pickled_tasks: FilePath) -> None:
    """Execute the tasks."""
    MPI, COMM = get_mpi()
    if COMM.rank == 0:
        with open(pickled_tasks, "rb") as pkl:
            tasks = pickle.load(pkl)
        mpi_jobs = split_tasks(tasks, COMM.size)
    else:
        mpi_jobs = None

    jobs: list[CNSJob] = COMM.scatter(mpi_jobs, root=0)

    results: list[FilePath] = []
    for job in jobs:
        job.run()
        # check if the job has an input file
        if hasattr(job, "input_file"):
            results.append(job.input_file)

    # COMM.Barrier()
    # results = MPI.COMM_WORLD.gather(results, root=0)


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
