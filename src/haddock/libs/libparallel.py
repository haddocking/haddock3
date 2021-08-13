"""Module in charge of parallelizing the execution of tasks"""
import logging
from multiprocessing import Process

from haddock.libs.libutil import parse_ncores


logger = logging.getLogger(__name__)


class Worker(Process):

    def __init__(self, tasks):
        super(Worker, self).__init__()
        self.tasks = tasks
        logger.info(f"Worker ready with {len(self.tasks)} tasks")

    def run(self):
        for task in self.tasks:
            task.run()
        logger.info(f"{self.name} executed")


class Scheduler:

    def __init__(self, tasks, ncores=None):
        """
        Schedule tasks to a defined number of processes.

        Parameters
        ----------
        tasks : list
            The list of tasks to execute. Tasks must be subclass of
            `multiprocessing.Process`.

        ncores : None or int
            The number of cores to use. If `None` is given uses the
            maximum number of CPUs allowed by
            `libs.libututil.parse_ncores` function.
        """

        self.num_tasks = len(tasks)
        self.num_processes = ncores  # first parses num_cores

        # Do not waste resources
        self.num_processes = min(self.num_processes, self.num_tasks)

        # step trick by @brianjimenez
        _n = self.num_processes
        job_list = [tasks[i::_n] for i in range(_n)]

        self.task_list = [Worker(jobs) for jobs in job_list]

        logger.info(f"{self.num_tasks} tasks ready.")

    @property
    def num_processes(self):
        return self._ncores

    @num_processes.setter
    def num_processes(self, n):
        self._ncores = parse_ncores(n)
        logger.info(f"Scheduler configurated for {self._ncores} cpu cores.")

    def run(self):

        try:
            for task in self.task_list:
                task.start()

            for task in self.task_list:
                task.join()

            logger.info(f"{self.num_tasks} tasks finished")

        except KeyboardInterrupt:
            # Q: why have a keyboard interrupt here?
            self.terminate()

    def terminate(self):

        logger.warning("Something went wrong")
        for task in self.task_list:
            task.terminate()

        logger.warning("The workers have stopped")
