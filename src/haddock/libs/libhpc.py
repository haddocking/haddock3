"""Module in charge of running tasks in HPC."""
from haddock import log
from pathlib import Path
import subprocess
import os
import shlex
import re
import time

STATE_REGEX = r"JobState=(\w*)"


class HPCWorker:
    """Defines the HPC Job."""

    def __init__(self, tasks, num):
        self.tasks = tasks
        log.debug(f"HPCWorker ready with {len(self.tasks)}")
        self.job_id = None
        self.job_num = num
        self.job_status = "unknown"
        self.job_fname = ""

    def run(self):
        """Execute the tasks."""
        job_file = self.prepare_job_file(self.tasks)
        cmd = f"sbatch {job_file}"
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
            # TODO: Maybe a regex here is overkill
            # https://regex101.com/r/M2vbAc/1
            status = re.findall(STATE_REGEX, out)[0]
            if status == "PENDING":
                self.job_status = "submitted"
            elif status == "RUNNING":
                self.job_status = "running"
            elif status == "SUSPENDED":
                self.job_status = "hold"
            elif status == "COMPLETING":
                self.job_status = "running"
            elif status == "COMPLETED":
                self.job_status = "finished"
            elif status == "FAILED":
                self.job_status = "failed"
        else:
            self.job_status = "finished"

        return self.job_status

    def prepare_job_file(self, job_list):
        job_name = "haddock3"
        queue = "haddock"
        moddir = job_list[0].modpath
        module = job_list[0].cns_folder
        run = job_list[0].config_path
        toppar = job_list[0].toppar
        module_name = job_list[0].modpath.name.split("_")[-1]
        self.job_fname = Path(moddir, f"{module_name}_{self.job_num}.job")
        out_fname = Path(moddir, f"{module_name}_{self.job_num}.out")
        err_fname = Path(moddir, f"{module_name}_{self.job_num}.err")

        header = f"#!/bin/sh{os.linesep}"
        header += f"#SBATCH -J {job_name}{os.linesep}"
        header += f"#SBATCH -p {queue}{os.linesep}"
        header += f"#SBATCH --nodes=1{os.linesep}"
        header += f"#SBATCH --tasks-per-node=1{os.linesep}"
        header += f"#SBATCH --output={out_fname}{os.linesep}"
        header += f"#SBATCH --error={err_fname}{os.linesep}"
        header += f"#SBATCH --workdir={moddir}{os.linesep}"

        envs = f"export MODDIR={moddir}{os.linesep}"
        envs += f"export MODULE={module}{os.linesep}"
        envs += f"export RUN={run}{os.linesep}"
        envs += f"export TOPPAR={toppar}{os.linesep}"

        body = envs

        body += f"cd {moddir}" + os.linesep
        for job in job_list:
            cmd = (
                f"{job.cns_exec} < {job.input_file} > {job.output_file}"
                f"{os.linesep}"
                )
            body += cmd

        with open(self.job_fname, "w") as job_fh:
            job_fh.write(header)
            job_fh.write(body)

        return self.job_fname

    def cancel(self):
        """Cancel the execution."""
        bypass_statuses = ["finished", "failed"]
        if self.update_status() not in bypass_statuses:
            log.info(f"Canceling {self.job_fname.name} - {self.job_id}")
            cmd = f"scancel {self.job_id}"
            _ = subprocess.run(shlex.split(cmd), capture_output=True)


class HPCScheduler:
    """Schedules tasks to run in HPC."""

    def __init__(self, task_list, queue_limit, concat):
        self.num_tasks = len(task_list)
        # FIXME: The defaults are hardcoded here
        if not concat:
            concat = 1
        if not queue_limit:
            queue_limit = 100
        # =======
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

        log.debug(f"{self.num_tasks} HPC tasks ready.")

    def run(self):
        """Run tasks in the Queue."""
        # split by maximum number of submission so we do it in batches
        adaptive_l = []
        batch = [
            self.worker_list[i:i + self.queue_limit]
            for i in range(0, len(self.worker_list), self.queue_limit)
            ]
        try:
            for batch_num, worker_list in enumerate(batch, start=1):
                log.info(f"> Running batch {batch_num}/{len(batch)}")
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

                    per = (
                        (completed_count + failed_count) / self.num_tasks
                        ) * 100
                    log.info(f">> {per:.0f}% done")
                    if completed_count + failed_count == len(worker_list):
                        completed = True
                        end = time.time()
                        elapsed = end - start
                        log.info(f">> Took {elapsed:.2f}s")
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
                log.info(f"> Batch {batch_num}/{len(batch)} done")

        except KeyboardInterrupt as err:
            self.terminate()
            raise err

    def terminate(self):
        """Terminate all jobs in the queue in a controlled way."""
        log.info("Terminate signal recieved, removing jobs from the queue...")
        for worker in self.worker_list:
            worker.cancel()

        log.info("The jobs in the queue were terminated in a controlled way")
