"""Module in charge of parallelizing the execution of tasks"""

import logging
from multiprocessing import Process, cpu_count

logger = logging.getLogger(__name__)


class Sailor(Process):
    """What is a captain without sailors?"""
    def __init__(self, tasks):
        super(Sailor, self).__init__()
        self.tasks = tasks
        logger.info(f'Sailor ready with {len(self.tasks)} tasks')

    def run(self):
        for task in self.tasks:
            task.run()
        logger.info(f'Sailor {self.name} says: Aye, Captain!')


class CaptainHaddock(object):
    """Aye, Captain"""
    def __init__(self, tasks, num_cores=0):
        try:
            self.num_processes = int(num_cores)
            if self.num_processes < 1:
                raise ValueError()
            if self.num_processes > cpu_count():
                logger.warning(f'Number of cores ({self.num_processes}) larger than available.')
                raise ValueError()
        except (ValueError, TypeError):
            logger.warning("Number of cores has not been specified or is incorrect. Using available cores.")
            self.num_processes = cpu_count()

        self.tasks = tasks
        self.num_tasks = len(tasks)
        # Do not waste resources
        if self.num_tasks < self.num_processes:
            self.num_processes = self.num_tasks
        logger.info(f'Captain commands {self.num_processes} sailors (cpu cores)')

        self.crew = []
        crew_tasks = [tasks[i::self.num_processes] for i in range(self.num_processes)]

        for i in range(self.num_processes):
            sailor = Sailor(crew_tasks[i])
            self.crew.append(sailor)

        logger.info(f'{self.num_tasks} tasks ready')

    def drink(self):
        """Drink up me 'earties, yo ho"""
        logger.info("Oh Captain, my Captain!")
        try:
            for sailor in self.crew:
                sailor.start()

            for sailor in self.crew:
                sailor.join()

            logger.info(f'{self.num_tasks} tasks finished')
        except KeyboardInterrupt:
            self.mutiny()

    def mutiny(self):
        """There is a mutiny in the ship, anyone interrupted?"""
        logger.warning("Mutiny in the ship")
        for sailor in self.crew:
            sailor.terminate()
        logger.warning("The crew has stopped working")
