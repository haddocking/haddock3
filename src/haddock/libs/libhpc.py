"""Module in charge of running tasks in HPC."""
import os
import re
import shlex
import subprocess
import time
from pathlib import Path

from haddock import log, modules_defaults_path
from haddock.gear.yaml2cfg import read_from_yaml_config


STATE_REGEX = r"JobState=(\w*)"

JOB_STATUS_DIC = {
    "PENDING": "submitted",
    "RUNNING": "running",
    "SUSPENDED": "hold",
    "COMPLETING": "running",
    "COMPLETED": "finished",
    "FAILED": "failed",
    }

# if you change these defaults, chage also the values in the
# modules/defaults.cfg file
_tmpcfg = read_from_yaml_config(modules_defaults_path)
HPCScheduler_CONCAT_DEFAULT = _tmpcfg["concat"]  # original value 1
HPCWorker_QUEUE_LIMIT_DEFAULT = _tmpcfg["queue_limit"]  # original value 100
HPCWorker_QUEUE_DEFAULT = _tmpcfg["queue"]  # original value ""
del _tmpcfg


class HPCWorker:
    """Defines the HPC Job."""

    def __init__(
            self,
            tasks,
            num,
            job_id=None,
            workfload_manager='slurm',
            queue=None
            ):
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

        self.moddir = Path(tasks[0].envvars['MODDIR'])
        self.toppar = tasks[0].envvars['TOPPAR']
        self.cns_folder = tasks[0].envvars['MODULE']
        module_name = \
            Path(tasks[0].envvars['MODDIR']).resolve().stem.split('_')[-1]
        self.job_fname = Path(self.moddir, f'{module_name}_{num}.job')
        self.workload_manager = workfload_manager
        self.queue = queue

    def prepare_job_file(self, queue_type='slurm'):
        """Prepare the job file for all the jobs in the task list."""
        job_file_contents = create_job_header_funcs[queue_type](
            job_name='haddock3',
            queue=self.queue,
            ncores=1,
            work_dir=self.moddir,
            stdout_path=self.job_fname.with_suffix('.out'),
            stderr_path=self.job_fname.with_suffix('.err'),
            )

        job_file_contents += create_CNS_export_envvars(
            MODDIR=self.moddir,
            MODULE=self.cns_folder,
            TOPPAR=self.toppar,
            )

        job_file_contents += f"cd {self.moddir}{os.linesep}"
        for job in self.tasks:
            cmd = (
                f"{job.cns_exec} < {job.input_file} > {job.output_file}"
                f"{os.linesep}"
                )
            job_file_contents += cmd

        self.job_fname.write_text(job_file_contents)

    def run(self):
        """Execute the tasks."""
        self.prepare_job_file(self.workload_manager)
        cmd = f"sbatch {self.job_fname}"
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        self.job_id = int(p.stdout.decode("utf-8").split()[-1])
        self.job_status = "submitted"

    def update_status(self):
        """Retrieve the status of this worker."""
        cmd = f"scontrol show jobid -dd {self.job_id}"
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        out = p.stdout.decode("utf-8")
        # err = p.stderr.decode('utf-8')
        if out:
            # https://regex101.com/r/M2vbAc/1
            status = re.findall(STATE_REGEX, out)[0]
            self.job_status = JOB_STATUS_DIC[status]
        else:
            self.job_status = "finished"

        return self.job_status

    def cancel(self, bypass_statuses=("finished", "failed")):
        """Cancel the execution."""
        if self.update_status() not in bypass_statuses:
            log.info(f"Canceling {self.job_fname.name} - {self.job_id}")
            cmd = f"scancel {self.job_id}"
            _ = subprocess.run(shlex.split(cmd), capture_output=True)


class HPCScheduler:
    """Schedules tasks to run in HPC."""

    def __init__(
            self,
            task_list,
            target_queue=HPCWorker_QUEUE_DEFAULT,
            queue_limit=HPCWorker_QUEUE_LIMIT_DEFAULT,
            concat=HPCScheduler_CONCAT_DEFAULT,
            ):
        self.num_tasks = len(task_list)
        self.queue_limit = queue_limit
        self.concat = concat

        # split tasks according to concat level
        if concat > 1:
            log.info(
                f"Concatenating, each .job will produce {concat} "
                "(or less) models"
                )
        job_list = [
            task_list[i:i + concat] for i in range(0, len(task_list), concat)
            ]

        self.worker_list = [
            HPCWorker(t, j) for j, t in enumerate(job_list, start=1)
            ]

        # set the queue
        #  (this is outside the comprehension for clarity)
        if target_queue:
            for worker in self.worker_list:
                worker.queue = target_queue

        log.debug(f"{self.num_tasks} HPC tasks ready.")

    def run(self):
        """Run tasks in the Queue."""
        # split by maximum number of submission so we do it in batches
        adaptive_l = []
        batch = [
            self.worker_list[i:i + self.queue_limit]
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
                completed = False
                while not completed:
                    for worker in worker_list:
                        worker.update_status()
                        if worker.job_status != "finished":
                            log.info(
                                f">> {worker.job_fname.name}"
                                f" {worker.job_status}"
                                )

                    completed_count = sum(
                        w.job_status == "finished" for w in worker_list
                        )

                    failed_count = sum(
                        w.job_status == "failed" for w in worker_list
                        )

                    if completed_count + failed_count == len(worker_list):
                        completed = True
                        end = time.time()
                        elapsed = end - start
                        adaptive_l.append(elapsed)
                    else:
                        if not adaptive_l:
                            # This is the first run, use pre-defined waits
                            if len(worker_list) < 10:
                                sleep_timer = 10
                            elif len(worker_list) < 50:
                                sleep_timer = 30
                            else:
                                sleep_timer = 60
                        else:
                            # We already know how long it took, use the average
                            sleep_timer = round(
                                sum(adaptive_l) / len(adaptive_l)
                                )
                        log.info(f">> Waiting... ({sleep_timer:.2f}s)")
                        time.sleep(sleep_timer)

                per = (float(batch_num) / float(total_batches)) * 100
                log.info(
                    f">> Batch {batch_num}/{total_batches} took "
                    f"{elapsed:.2f}s to finish, {per:.2f}% complete")

        except KeyboardInterrupt as err:
            self.terminate()
            raise err

    def terminate(self):
        """Terminate all jobs in the queue in a controlled way."""
        log.info("Terminate signal recieved, removing jobs from the queue...")
        for worker in self.worker_list:
            worker.cancel()

        log.info("The jobs in the queue were terminated in a controlled way")


def create_slurm_header(
        job_name='haddock3_slurm_job',
        work_dir='.',
        stdout_path='haddock3_job.out',
        stderr_path='haddock3_job.err',
        queue=None,
        ncores=48,
        ):
    """
    Create HADDOCK3 Slurm Batch job file.

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
    header += f"#SBATCH --workdir={work_dir}{os.linesep}"
    return header


def create_torque_header(
        job_name='haddock3_torque_job',
        work_dir='.',
        stdout_path='haddock3_job.out',
        stderr_path='haddock3_job.err',
        queue=None,
        ncores=48,
        ):
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


def create_CNS_export_envvars(**envvars):
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
        f'export {key.upper()}={value}'
        for key, value in envvars.items()
        )

    return exports + os.linesep + os.linesep


# the different job submission queues
create_job_header_funcs = {
    'torque': create_torque_header,
    'slurm': create_slurm_header,
    }
