"""Module in charge of parallelizing the execution of tasks"""

import logging
from multiprocessing import Process, cpu_count

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


class Scheduler():

    def __init__(self, tasks, num_cores=0):
        try:
            self.num_processes = int(num_cores)

            if self.num_processes < 1:
                raise ValueError()

            if self.num_processes > cpu_count():
                logger.warning(f"Number of cores ({self.num_processes}) larger"
                               " than available.")
                raise ValueError()

        except (ValueError, TypeError):
            logger.warning("Number of cores has not been specified or is"
                           " incorrect. Using available cores.")
            self.num_processes = cpu_count()

        self.tasks = tasks
        self.num_tasks = len(tasks)
        # Do not waste resources
        if self.num_tasks < self.num_processes:
            self.num_processes = self.num_tasks
        logger.info(f"Running with {self.num_processes} cpu cores")

        self.task_list = []
        job_list = []
        for i in range(self.num_processes):
            job_list.append(tasks[i::self.num_processes])

        for i in range(self.num_processes):
            task = Worker(job_list[i])
            self.task_list.append(task)

        logger.info(f"{self.num_tasks} tasks ready")

    def execute(self):

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
