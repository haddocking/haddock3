"""Module in charge of parallelizing the execution of tasks"""
from multiprocessing import Process

from haddock import log
from haddock.libs.libutil import parse_ncores


class Worker(Process):

    def __init__(self, tasks):
        super(Worker, self).__init__()
        self.tasks = tasks
        log.info(f"Worker ready with {len(self.tasks)} tasks")

    def run(self):
        for task in self.tasks:
            task.run()
        log.info(f"{self.name} executed")


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

        log.info(f"{self.num_tasks} tasks ready.")

    @property
    def num_processes(self):
        return self._ncores

    @num_processes.setter
    def num_processes(self, n):
        self._ncores = parse_ncores(n)
        log.info(f"Scheduler configurated for {self._ncores} cpu cores.")

    def run(self):

        try:
            for task in self.task_list:
                task.start()

            for task in self.task_list:
                task.join()

            log.info(f"{self.num_tasks} tasks finished")

        except KeyboardInterrupt as err:
            # Q: why have a keyboard interrupt here?
            # A: To have a controlled break if the user Ctrl+c during CNS run
            self.terminate()
            # this raises sends the error to libs.libworkflow.Step
            # if Scheduler is used independently the error will propagate to
            # whichever has to catch it
            raise err

    def terminate(self):

        for task in self.task_list:
            task.terminate()

        log.info("The workers terminated in a controlled way")
