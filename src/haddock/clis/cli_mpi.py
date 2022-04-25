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


try:
    from mpi4py import MPI
except ImportError as e:
    _msg = (
        f"{e} - To run this cli you must have mpi4py and "
        "OpenMPI installed in the system")
    sys.exit(_msg)


COMM = MPI.COMM_WORLD


def split_tasks(task_l, n):
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


def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(ap, main)


# ========================================================================#


def main(pickled_tasks):
    """Execute the tasks."""
    if COMM.rank == 0:
        with open(pickled_tasks, "rb") as pkl:
            tasks = pickle.load(pkl)
        jobs = split_tasks(tasks, COMM.size)
    else:
        jobs = None

    jobs = COMM.scatter(jobs, root=0)

    results = []
    for job in jobs:
        job.run()
        results.append(job.input_file)

    # COMM.Barrier()
    # results = MPI.COMM_WORLD.gather(results, root=0)


if __name__ == "__main__":
    sys.exit(maincli())
