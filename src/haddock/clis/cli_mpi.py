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


try:
    from mpi4py import MPI
except ImportError as e:
    _msg = (
        f"{e} - To run this cli you must have mpi4py and "
        "OpenMPI installed in the system"
    )
    sys.exit(_msg)


COMM = MPI.COMM_WORLD


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
        results.append(job.input_file)

    # COMM.Barrier()
    # results = MPI.COMM_WORLD.gather(results, root=0)


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
