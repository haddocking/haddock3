"""Module in charge of parallelizing the execution of tasks."""

import math
from multiprocessing import Process, Queue

from haddock import log
from haddock.core.typing import (
    AnyT,
    FilePath,
    Generator,
    Optional,
    Sequence,
    SupportsRunT,
    Union,
)
from haddock.libs.libutil import parse_ncores


def split_tasks(lst: Sequence[AnyT], n: int) -> Generator[Sequence[AnyT], None, None]:
    """Split tasks into N-sized chunks."""
    n = math.ceil(len(lst) / n)
    for j in range(0, len(lst), n):
        chunk = lst[j : n + j]
        yield chunk


def get_index_list(nmodels, ncores):
    """
    Optimal distribution of models among cores

    Parameters
    ----------
    nmodels : int
        Number of models to be distributed.

    ncores : int
        Number of cores to be used.

    Returns
    -------
    index_list : list
        List of model indexes to be used for the parallel scanning.
    """
    if nmodels < 1:
        raise ValueError(f"nmodels ({nmodels})) must be greater than 0")
    if ncores < 1:
        raise ValueError(f"ncores ({ncores}) must be greater than 0")
    spc = nmodels // ncores
    # now the remainder
    rem = nmodels % ncores
    # now the list of indexes to be used for the SCAN calculation
    index_list = [0]
    for core in range(ncores):
        if core < rem:
            index_list.append(index_list[-1] + spc + 1)
        else:
            index_list.append(index_list[-1] + spc)
    return index_list


class GenericTask:
    """Generic task to be executed."""

    def __init__(self, function, *args, **kwargs):
        if not callable(function):
            raise TypeError("The 'function' argument must be callable")
        self.function = function
        self.args = args
        self.kwargs = kwargs

    def run(self):
        return self.function(*self.args, **self.kwargs)


class Worker(Process):
    """Work on tasks."""

    def __init__(self, tasks: Sequence[SupportsRunT], results: Queue) -> None:
        super(Worker, self).__init__()
        self.tasks = tasks
        self.result_queue = results
        log.debug(f"Worker ready with {len(self.tasks)} tasks")

    def run(self) -> None:
        """Execute tasks."""
        results = []
        for task in self.tasks:
            r = None
            try:
                r = task.run()
            except Exception as e:
                log.warning(f"Exception in task execution: {e}")

            results.append(r)

        # Put results into the queue
        self.result_queue.put(results)

        # Signal completion by putting a unique identifier into the queue
        self.result_queue.put(f"{self.name}_done")

        # log.debug(f"{self.name} executed")


class Scheduler:
    """Schedules tasks to run in multiprocessing."""

    def __init__(
        self,
        tasks: list[SupportsRunT],
        ncores: Optional[int] = None,
        max_cpus: bool = False,
    ) -> None:
        """
        Schedule tasks to a defined number of processes.

        Parameters
        ----------
        tasks : list
            The list of tasks to execute. Tasks must have method `run()`.

        ncores : None or int
            The number of cores to use. If `None` is given uses the
            maximum number of CPUs allowed by
            `libs.libututil.parse_ncores` function.
        """
        self.max_cpus = max_cpus
        self.num_tasks = len(tasks)
        self.num_processes = ncores  # first parses num_cores
        self.queue: Queue = Queue()
        self.results: list = []

        # Sort the tasks by input_file name and its length, so we know that 2 comes before 10
        ### Q? Whys is this necessary?
        # Only CNSJobs can be sorted like this
        if all(hasattr(t, "input_file") for t in tasks):
            task_name_dic: dict[int, tuple[FilePath, int]] = {}
            for i, t in enumerate(tasks):
                task_name_dic[i] = (t.input_file, len(str(t.input_file)))  # type: ignore

            sorted_task_list: list[SupportsRunT] = []
            for e in sorted(task_name_dic.items(), key=lambda x: (x[0], x[1])):
                idx = e[0]
                sorted_task_list.append(tasks[idx])
        else:
            sorted_task_list = tasks

        job_list = split_tasks(sorted_task_list, self.num_processes)
        self.worker_list = [Worker(jobs, self.queue) for jobs in job_list]

        log.info(f"Using {self.num_processes} cores")
        log.debug(f"{self.num_tasks} tasks ready.")

    @property
    def num_processes(self) -> int:
        """Number of processors to use."""  # noqa: D401
        return self._ncores

    @num_processes.setter
    def num_processes(self, n: Union[str, int, None]) -> None:
        self._ncores = parse_ncores(
            n,
            njobs=self.num_tasks,
            max_cpus=self.max_cpus,
        )
        log.debug(f"Scheduler configured for {self._ncores} cpu cores.")

    def run(self) -> None:
        """Run tasks in parallel."""

        try:
            for w in self.worker_list:
                w.start()

            # Collect results until all workers have signaled completion
            all_results = []
            num_workers = len(self.worker_list)
            completed_workers = 0

            while completed_workers < num_workers:
                result = self.queue.get()
                if isinstance(result, str) and result.endswith("_done"):
                    completed_workers += 1
                else:
                    all_results.append(result)

            for w in self.worker_list:
                w.join()

            self.results = [item for sublist in all_results for item in sublist]

            log.info(f"{self.num_tasks} tasks finished")

        except KeyboardInterrupt as err:
            # Q: why have a keyboard interrupt here?
            # A: To have a controlled break if the user Ctrl+c during CNS run
            self.terminate()
            # this raises sends the error to libs.libworkflow.Step
            # if Scheduler is used independently the error will propagate to
            # whichever has to catch it
            raise err

    def terminate(self) -> None:
        """Terminate tasks in a controlled way."""
        for worker in self.worker_list:
            worker.terminate()

        log.info("The workers terminated in a controlled way")
