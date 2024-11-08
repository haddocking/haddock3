"""Module in charge of running tasks in HPC."""

import os
import re
import shlex
import subprocess
import time
from pathlib import Path

from haddock import log, modules_defaults_path
from haddock.core.typing import Any, Container, FilePath, Optional
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs.libsubprocess import CNSJob


STATE_REGEX = r"JobState=(\w*)"

JOB_STATUS_DIC = {
    "PENDING": "submitted",
    "RUNNING": "running",
    "SUSPENDED": "hold",
    "COMPLETING": "running",
    "COMPLETED": "finished",
    "FAILED": "failed",
    "TIMEOUT": "timed-out",
}

TERMINATED_STATUS = (
    "finished",
    "failed",
    "timed-out",
)

# if you change these defaults, change also the values in the
# modules/defaults.cfg file
_tmpcfg = read_from_yaml_config(modules_defaults_path)
HPCScheduler_CONCAT_DEFAULT: int = _tmpcfg["concat"]  # original value 1
HPCWorker_QUEUE_LIMIT_DEFAULT: int = _tmpcfg[
    "queue_limit"
]  # original value 100 # noqa: E501
HPCWorker_QUEUE_DEFAULT: str = _tmpcfg["queue"]  # original value ""
del _tmpcfg


class HPCWorker:
    """Defines the HPC Job."""

    def __init__(
        self,
        tasks: list[CNSJob],
        num: int,
        job_id: Optional[int] = None,
        workfload_manager: str = "slurm",
        queue: Optional[str] = None,
    ) -> None:
        """
        Define the HPC job.

        Parameters
        ----------
        tasks : list of libs.libcns.CNSJob objects

        num : int
            The number of the worker.
        """
        self.tasks = tasks
        log.debug(f"HPCWorker ready with {len(self.tasks)}")
        self.job_num = num
        self.job_id = job_id
        self.job_status = "unknown"

        self.moddir = Path(tasks[0].envvars["MODDIR"])
        self.toppar = tasks[0].envvars["TOPPAR"]
        self.cns_folder = tasks[0].envvars["MODULE"]
        module_name = Path(tasks[0].envvars["MODDIR"]).resolve().stem.split("_")[-1]
        self.job_fname = Path(self.moddir, f"{module_name}_{num}.job")
        self.workload_manager = workfload_manager
        self.queue = queue

    def prepare_job_file(self, queue_type: str = "slurm") -> None:
        """Prepare the job file for all the jobs in the task list."""
        job_file_contents = create_job_header_funcs[queue_type](
            job_name="haddock3",
            queue=self.queue,
            ncores=1,
            work_dir=self.moddir,
            stdout_path=self.job_fname.with_suffix(".out"),
            stderr_path=self.job_fname.with_suffix(".err"),
        )

        job_file_contents += create_CNS_export_envvars(
            MODDIR=self.moddir,
            MODULE=self.cns_folder,
            TOPPAR=self.toppar,
        )

        job_file_contents += f"cd {self.moddir}{os.linesep}"
        for job in self.tasks:
            cmd = (
                f"{job.cns_exec} < {job.input_file} > {job.output_file}" f"{os.linesep}"
            )
            job_file_contents += cmd

        self.job_fname.write_text(job_file_contents)

    def run(self) -> None:
        """Execute the tasks."""
        self.prepare_job_file(queue_type=self.workload_manager)
        cmd = f"sbatch {self.job_fname}"
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        self.job_id = int(p.stdout.decode("utf-8").split()[-1])
        self.job_status = "submitted"

    def update_status(self) -> str:
        """Retrieve the status of this worker."""
        cmd = f"scontrol show jobid -dd {self.job_id}"
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        out = p.stdout.decode("utf-8")
        # err = p.stderr.decode('utf-8')
        if out:
            status = extract_slurm_status(out)
            self.job_status = JOB_STATUS_DIC[status]
        else:
            self.job_status = "finished"

        return self.job_status

    def cancel(self, bypass_statuses: Container[str] = ("finished", "failed")) -> None:
        """Cancel the execution."""
        if self.update_status() not in bypass_statuses:
            log.info(f"Canceling {self.job_fname.name} - {self.job_id}")
            cmd = f"scancel {self.job_id}"
            _ = subprocess.run(shlex.split(cmd), capture_output=True)


class HPCScheduler:
    """Schedules tasks to run in HPC."""

    def __init__(
        self,
        task_list: list[CNSJob],
        target_queue: str = HPCWorker_QUEUE_DEFAULT,
        queue_limit: int = HPCWorker_QUEUE_LIMIT_DEFAULT,
        concat: int = HPCScheduler_CONCAT_DEFAULT,
    ) -> None:
        self.num_tasks = len(task_list)
        self.queue_limit = queue_limit
        self.concat = concat

        # split tasks according to concat level
        if concat > 1:
            log.info(
                f"Concatenating, each .job will produce {concat} " "(or less) models"
            )
        job_list = [task_list[i : i + concat] for i in range(0, len(task_list), concat)]

        self.worker_list = [HPCWorker(t, j) for j, t in enumerate(job_list, start=1)]

        # set the queue
        #  (this is outside the comprehension for clarity)
        if target_queue:
            for worker in self.worker_list:
                worker.queue = target_queue

        log.debug(f"{self.num_tasks} HPC tasks ready.")

    def run(self) -> None:
        """Run tasks in the Queue."""
        # split by maximum number of submission so we do it in batches
        batch = [
            self.worker_list[i : i + self.queue_limit]
            for i in range(0, len(self.worker_list), self.queue_limit)
        ]
        total_batches = len(batch)
        try:
            for batch_num, worker_list in enumerate(batch, start=1):
                log.info(f"> Running batch {batch_num}/{total_batches}")
                start = time.time()
                for worker in worker_list:
                    worker.run()

                # check if those finished
                completed: bool = False
                elapsed: float = 0.0
                while not completed:
                    # Initiate count of terminated jobs
                    terminated_count: int = 0
                    # Loop over workers
                    for worker in worker_list:
                        worker.update_status()
                        # Log status if not finished
                        if worker.job_status != "finished":
                            log.info(
                                f">> {worker.job_fname.name}" f" {worker.job_status}"
                            )
                        # Increment number of terminated works
                        if worker.job_status in TERMINATED_STATUS:
                            terminated_count += 1

                    # Check if all terminated
                    if terminated_count == len(worker_list):
                        # Set while loop condition
                        completed = True
                        end = time.time()
                        elapsed = end - start
                    else:
                        # use pre-defined waits
                        if len(worker_list) < 10:
                            sleep_timer = 10
                        elif len(worker_list) < 50:
                            sleep_timer = 30
                        else:
                            sleep_timer = 60
                        log.info(f">> Waiting... ({sleep_timer:.2f}s)")
                        time.sleep(sleep_timer)

                per = (float(batch_num) / float(total_batches)) * 100
                log.info(
                    f">> Batch {batch_num}/{total_batches} took "
                    f"{elapsed:.2f}s to finish, {per:.2f}% complete"
                )

        except KeyboardInterrupt as err:
            self.terminate()
            raise err

    def terminate(self) -> None:
        """Terminate all jobs in the queue in a controlled way."""
        log.info("Terminate signal received, removing jobs from the queue...")
        for worker in self.worker_list:
            worker.cancel()

        log.info("The jobs in the queue were terminated in a controlled way")


def create_slurm_header(
    job_name: FilePath = "haddock3_slurm_job",
    work_dir: FilePath = ".",
    stdout_path: FilePath = "haddock3_job.out",
    stderr_path: FilePath = "haddock3_job.err",
    queue: Optional[str] = None,
    ncores: int = 48,
) -> str:
    """
    Create HADDOCK3 Slurm Batch job file.

    Parameters
    ----------
    job_name : str
        The name of the job.

    work_dir : pathlib.Path
        The working dir of the example. That is, the directory where
        `input`, `jobs`, and `logs` reside. Injected in `create_job_header`.

    time : int
        Time in minutes before job reach TIMEOUT status.

    **job_params
        According to `job_setup`.

    Return
    ------
    str
        Slurm-based job file for HADDOCK3.
    """
    header = f"#!/usr/bin/env bash{os.linesep}"
    header += f"#SBATCH -J {job_name}{os.linesep}"
    if queue:
        header += f"#SBATCH -p {queue}{os.linesep}"
    header += f"#SBATCH --nodes=1{os.linesep}"
    header += f"#SBATCH --tasks-per-node={str(ncores)}{os.linesep}"
    header += f"#SBATCH --output={stdout_path}{os.linesep}"
    header += f"#SBATCH --error={stderr_path}{os.linesep}"
    # commenting the workdir option (not supported by all versions of slurm)
    # header += f"#SBATCH --workdir={work_dir}{os.linesep}"
    return header


def create_torque_header(
    job_name: FilePath = "haddock3_slurm_job",
    work_dir: FilePath = ".",
    stdout_path: FilePath = "haddock3_job.out",
    stderr_path: FilePath = "haddock3_job.err",
    queue: Optional[str] = None,
    ncores: int = 48,
) -> str:
    """
    Create HADDOCK3 Alcazar job file.

    Parameters
    ----------
    job_name : str
        The name of the job.

    work_dir : pathlib.Path
        The working dir of the example. That is, the directory where
        `input`, `jobs`, and `logs` reside. Injected in `create_job_header`.

    **job_params
        According to `job_setup`.

    Return
    ------
    str
        Torque-based job file for HADDOCK3 benchmarking.
    """
    header = f"#!/usr/bin/env tcsh{os.linesep}"
    header += f"#PBS -N {job_name}{os.linesep}"
    if queue:
        header += f"#PBS -q {queue}{os.linesep}"
    header += f"#PBS -l nodes=1:ppn={str(ncores)}{os.linesep}"
    header += f"#PBS -S /bin/tcsh{os.linesep}"
    header += f"#PBS -o {stdout_path}{os.linesep}"
    header += f"#PBS -e {stderr_path}{os.linesep}"
    header += f"#PBS -wd {work_dir}{os.linesep}"
    return header


def to_torque_time(time: int) -> str:
    """Convert time in minutes to the form hh:mm:ss.

    Parameters
    ----------
    time : int
        Time in minutes.

    Return
    ------
    hh_mm_ss : str
        Time in the form for HH:MM:SS
    """
    hours = time // 60
    remain_mins = time - (hours * 60)
    # Convert to hh:mm:ss string
    hh_mm_ss_l = [hours, remain_mins, 0]
    # Make sure hours contain at least 2 characters
    hh_mm_ss = "{0:02d}:{1:02d}:{2:02d}".format(*hh_mm_ss_l)
    return hh_mm_ss


def extract_slurm_status(slurm_out: str) -> str:
    """Extract job status from slurm scontrol stdout.

    Parameters
    ----------
    slurm_out : str
        StdOut of `scontrol show jobid -dd {job_id}` command.
    Return
    ------
    status : str
        Status of the slurm job.
        May also return `error`, when job do not exists.
    """
    try:
        # https://regex101.com/r/M2vbAc/1
        status = re.findall(STATE_REGEX, slurm_out)[0]
    except IndexError:
        status = "FAILED"
    return status


def create_CNS_export_envvars(**envvars: Any) -> str:
    """Create a string exporting envvars needed for CNS.

    Parameters
    ----------
    envvars : dict
        A dictionary containing envvariables where keys are var names
        and values are the values.

    Returns
    -------
    str
        In the form of:
        export VAR1=VALUE1
        export VAR2=VALUE2
        export VAR3=VALUE3

    """
    exports = os.linesep.join(
        f"export {key.upper()}={value}" for key, value in envvars.items()
    )

    return exports + os.linesep + os.linesep


# the different job submission queues
create_job_header_funcs = {
    "torque": create_torque_header,
    "slurm": create_slurm_header,
}
