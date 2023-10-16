r"""
HADDOCK3 benchmark submission daemon.

For more information read our benchmark tutorial at `docs/benchmark.tut`
in HADDOCK3 repository site: https://github.com/haddocking/haddock3::

   (_) L|J
   (")  |
   /_\--|
 _/\ /  |
   _W_  |

Usage::

    haddock3-dmn -h
    haddock3-dmn <benchmark folder>  --job-limit <num> [OPTIONS]
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

from haddock.core.typing import ArgumentParser, Callable, Namespace, Optional


workload_manager_launch = {
    "slurm": "sbatch",
    "torque": "qsub",
}
"""options for the different job queue systems supported"""


# prepares client arguments
ap = argparse.ArgumentParser(
    prog="HADDOCK3 benchmark submission daemon.",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

ap.add_argument(
    "benchmark_path",
    help="Path to the benchmark folder as prepared by `haddock3-bm` interface.",
    type=Path,
)

ap.add_argument(
    "--job-limit",
    dest="job_limit",
    help="How many jobs should run at the same time. Default: 10",
    default=10,
    type=int,
)

ap.add_argument(
    "--job-sys",
    dest="manager",
    help="The system where the jobs will be run. Default `slurm`.",
    choices=tuple(workload_manager_launch.keys()),
    default="slurm",
)

ap.add_argument(
    "--restart",
    help="Restart the RUNNING jobs. DONE jobs won't be touched.",
    action="store_true",
)

ap.add_argument(
    "--sort-first",
    dest="sort_first",
    help=(
        "Sort jobs by size in ascending order. If not given jobs are order by "
        "size in descending order: the biggest first."
    ),
    action="store_true",
)


def _ap() -> ArgumentParser:
    return ap


class Job:
    """
    Job task.

    Controls the status of each job.

    Parameters
    ----------
    job_f : pathlib.Path
        The path to the job file.

    launch_command : str
        The command to launch the job. For example `sbatch`.
    """

    def __init__(self, job_f: Path, launch_command: str) -> None:
        self.job_filename = job_f

        self.launch_command = launch_command

        # previous jog path
        job_stem = job_f.stem
        self.job_run_folder = Path(job_f.parents[1], f"run-{job_stem}")

        self.check_done = Path(self.job_run_folder, "DONE")
        self.check_running = Path(self.job_run_folder, "RUNNING")
        self.check_available = Path(self.job_run_folder, "AVAILABLE")
        self.check_fail = Path(self.job_run_folder, "FAIL")
        self.status = None

        self.status_files = [
            self.check_done,
            self.check_running,
            self.check_fail,
            self.check_available,
        ]

    def get_status(self) -> Optional[str]:
        """
        Get job status.

        The job status depends on the present of files: `AVAILABLE`,
        `RUNNING`, `DONE`, and `FAIL` created by `haddock-bm` jobs.

        Status is assigned to `self.status` and returned.
        """
        for _file in self.status_files:
            if _file.exists():
                self.status = _file.stem  # type: ignore
                break
        return self.status

    def submit(self) -> None:
        """
        Submit job.

        Run command `$launch_command $job_filename`.
        """
        subprocess.run(cmds := [self.launch_command, str(self.job_filename)])
        print("Job sent: ", cmds)

    def restart(self) -> None:
        """
        Restart the status of the job to `AVAILABLE`.

        Does this by removing all status files and creating the file
        `AVAILABLE`.
        """
        self.check_done.unlink(missing_ok=True)
        self.check_running.unlink(missing_ok=True)
        self.check_fail.unlink(missing_ok=True)
        self.check_available.touch(exist_ok=True)


def get_current_jobs(grep: str = "BM5") -> int:
    """
    Get current number of jobs for which job-name has the `grep` word.

    List of jobs is retrieve using the command `qstat`.

    Parameters
    ----------
    grep : str
        The string to search job-names for.

    Return
    ------
    int
        The number of jobs with the word `grep` in their name.
    """
    concurrent_cmd = "qstat -a | awk '{print $4}' | grep " + grep + " | wc -l"
    p = subprocess.Popen(
        concurrent_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output = p.communicate()[0]
    njobs = int(output.decode("utf-8"))
    print(f"Found {njobs} in the queue.")
    return njobs


def calc_size(job_path: Path) -> int:
    """
    Calculate the size of the job.

    Expects the job file to be in a folder structure as defined by
    `haddock3-bm`.

    The size of the jobs is defined by the number of carbon-alpha lines
    in the `target.pdb` file.
    """
    target_pdb = Path(job_path.parents[1], "input", "target.pdb")
    lines = target_pdb.read_text().split(os.linesep)
    size = sum(1 for line in lines if line[11:16].strip() == "CA")
    return size


def filter_by_status(job_list: list[Job], status: str = "AVAILABLE") -> list[Job]:
    """
    Filter jobs by their status.

    Only jobs with `status` are accepted.

    Parameters
    ----------
    job_list : list of Job objects.

    Returns
    -------
    list
        The list with the `Job`s with matching `status`.
    """
    jobs = [j for j in job_list if j.get_status() == status]
    print(f"Number of {status} jobs: {len(jobs)}.")
    return jobs


# command-line client helper functions
# load_args, cli, maincli
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


def main(
    benchmark_path: Path,
    job_limit: int = 10,
    manager: str = "slurm",
    restart: bool = False,
    sort_first: bool = False,
) -> None:
    """
    Execute the benchmark daemon.

    The parameters defined here are the same as defined in the client
    arguments.

    This is the main function of the client. If you want to run the
    daemon withOUT using the command line and instead importing its
    functionalities and setting it up from another pythong script, you
    should import this function.

    >>> from haddock.clis.cli_dmn import main

    Parameters
    ----------
    benchmark_path : pathlib.Path
        The path of the benchmark folder as created by the `haddock3-bm`
        interface.

    job_limit : int
        The max number of jobs to send to the queue.

    manager : str
        A key to the `workload_manager_launch` dictionary. Selects the queue
        management system.

    restart : bool
        Whether to restart the `RUNNING` jobs that might have been halted
        in previous daemon runs. Defaults to False.

    sort_first : bool
        Whether to sort jobs by their size in ascending manner. That is,
        the sorted jobs first. Defaults to False, the longer first.
    """
    # lists of all the job files in the benchmark_path folder
    job_list = list(benchmark_path.glob("*/jobs/*.job"))

    # breaks if no jobs are found
    if not job_list:
        sys.exit("+ ERROR! No jobs found in folder: {str(benchmark_path)!r}")

    # sorts the job list
    job_list.sort(key=calc_size, reverse=not (sort_first))  # noqa: E275

    # create the job objects according to the queue managing systme
    _jobsys = workload_manager_launch[manager]
    jobs = [Job(j, _jobsys) for j in job_list]

    # restart previously (halted) `RUNNING` jobs - if selected.
    if restart:
        running_jobs = filter_by_status(jobs, status="RUNNING")
        for _job in running_jobs:
            _job.restart()

    # Lists the available jobs (those with status `AVAILABLE`)
    available_jobs = filter_by_status(jobs)

    # runs the daemon loop, only if there are available jobs :-)
    while available_jobs:
        # get the number of available queue slots according to the
        # job limit parameter
        #
        # 0 is added to avoid going to negative values in case a job
        # had been manually submitted.
        empty_slots = max(0, job_limit - get_current_jobs())
        print("empty slots: ", empty_slots)

        # send jobs if there are empty slots
        for i in range(empty_slots):
            available_jobs[i].submit()
            time.sleep(5)

        # chill before repeating the process.
        print("chilling for 120 seconds...")
        time.sleep(120)

        # refreshes the available_jobs list
        available_jobs = filter_by_status(jobs)

    # done
    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
