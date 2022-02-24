import argparse
import pickle
import sys

try:
    from mpi4py import MPI
except ImportError as e:
    sys.exit(e)


COMM = MPI.COMM_WORLD


def split_tasks(task_l, n):
    """Split tasks into equal lenght chunks."""
    return [task_l[_i::n] for _i in range(n)]


# @joaomcteixeira factory pattern loader stuff #=============================#
ap = argparse.ArgumentParser()

ap.add_argument(
    "pickled_tasks",
    help="The input pickled tasks path",
)


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
