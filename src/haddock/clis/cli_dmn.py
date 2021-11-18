r"""
HADDOCK3 benchmark submission daemon.

   (_) L|J
   (")  |
   /_\--|
 _/\ /  |
   _W_  |

Usage:
    haddock3-dmn <benchmark folder>  --job-limit <num>
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


job_system_launch = {
    'slurm': 'sbatch',
    'torque': 'qsub',
    }


ap = argparse.ArgumentParser(
    prog='HADDOCK3 benchmark submission daemon.',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    "benchmark_path",
    help='Path to the benchmark folder as prepared by `haddock3-bm` interface.',
    type=Path,
    )

ap.add_argument(
    "--job-limit",
    dest='job_limit',
    help="How many jobs should run at the same time. Default: 10",
    default=10,
    type=int,
    )

ap.add_argument(
    '--job-sys',
    dest='job_sys',
    help='The system where the jobs will be run. Default `slurm`.',
    choices=tuple(job_system_launch.keys()),
    default='slurm',
    )


class Job:
    """Job task."""

    def __init__(self, job_f, job_system):
        self.job_filename = job_f

        self.job_system = job_system

        # previous jog path
        job_stem = job_f.stem
        self.job_run_folder = Path(job_f.parents[1], f'run-{job_stem}')

        self.check_done = Path(self.job_run_folder, 'DONE')
        self.check_running = Path(self.job_run_folder, 'RUNNING')
        self.check_available = Path(self.job_run_folder, 'AVAILABLE')
        self.check_fail = Path(self.job_run_folder, 'FAIL')
        self.status = None

        self.status_files = [
            self.check_done,
            self.check_running,
            self.check_fail,
            self.check_available
            ]

    def get_status(self):
        """
        Get job status.

        If job has no status assigned in the disk, assigns 'AVAILABLE'.

        Status is assigned to `self.status` and returned.
        """
        for _file in self.status_files:
            if _file.exists():
                self.status = _file.stem
                break
        return self.status

    def submit(self):
        """Submit job."""
        subprocess.run([self.job_system, self.job_filename])
        print('Job sent: ', [self.job_system, self.job_filename])


def get_current_jobs(grep='BM5'):
    """
    Get current jobs for which job-name has the `grep` word.

    Parameters
    ----------
    grep : str
        The string to search job-names for.

    Return
    ------
    int
        The number of jobs with the word `grep` in their name.
    """
    concurrent_cmd = "qstat -a | awk '{print $4}' | grep ' + grep + ' | wc -l"
    p = subprocess.Popen(
        concurrent_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        )
    output = p.communicate()[0]
    return int(output.decode('utf-8'))


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


def calc_size(job_path):
    """
    Calculate the size of the job.

    Expects the job file to be in a folder structure as defined by
    `haddock3-bm`.
    """
    target_pdb = Path(job_path.parents[1], 'input', 'target.pdb')
    lines = target_pdb.read_text().split(os.linesep)
    size = sum(1 for line in lines if line[11:16].strip() == 'CA')
    return size


def filter_available(job_list, status='AVAILABLE'):
    """
    Filter jobs by their availability.

    Only jobs with `status` are accepted.

    Returns
    -------
    list
        The list with the `Job`s with matching `status`.
    """
    jobs = [j for j in job_list if j.get_status() == status]
    print(f'Number of {status} jobs: {len(jobs)}.')
    return jobs


def main(benchmark_path, job_limit=10, job_sys='slurm'):
    """."""
    job_list = list(benchmark_path.glob('*/jobs/*.job'))

    if not job_list:
        sys.exit('+ ERROR! No jobs found in folder: {str(benchmark_path)!r}')

    job_list.sort(key=calc_size, reverse=True)

    _jobsys = job_system_launch[job_sys]
    jobs = [Job(j, _jobsys) for j in job_list]
    available_jobs = filter_available(jobs)

    while available_jobs:

        # 0 is added to avoid going to negative values in case a job
        # had been manually submitted.
        empty_spots = max(0, job_limit - get_current_jobs())
        print('empty spots: ', empty_spots)

        for i in range(empty_spots):
            available_jobs[i].submit()
            time.sleep(5)

        # chill
        print('chilling for 120 seconds...')
        time.sleep(120)

        available_jobs = filter_available(jobs)

    return


if __name__ == '__main__':
    sys.exit(maincli())
