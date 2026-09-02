"""Module in charge of parallelizing the execution of tasks."""

import math
import queue
import threading
from concurrent.futures import Future, ThreadPoolExecutor
from multiprocessing import Process, Queue

from haddock import log
from haddock.core.exceptions import (
    CachedCNSFailure,
    ConfigurationError,
)
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
from haddock.libs.libseamless import (
    precompute_invariant_checksums_for_jobs,
    write_cns_dependencies,
)


_CACHE_WRITER_POLL_SECONDS = 0.05
_CACHE_WRITER_CHECKSUM_THREADS = 8


class CacheRecordWriter:
    """Write CNS cache records from one thread in the scheduler's process.

    Jobs are considered safe to checksum only after their worker has reported
    completion.  The writer still polls all outstanding output paths, which
    lets it record completed jobs during a long scheduler run without guessing
    at filesystem visibility latency.
    """

    def __init__(
        self,
        tasks: Sequence[SupportsRunT],
        completion_queue: Queue,
        scheduler_shutdown: threading.Event,
        checksum_workers: int = _CACHE_WRITER_CHECKSUM_THREADS,
    ) -> None:
        self.tasks = {
            task.cache_writer_id: task
            for task in tasks
            if getattr(task, "cache_context", None) is not None
        }
        self.completion_queue = completion_queue
        self.outstanding = set(self.tasks)
        self.completed: set[int] = set()
        self.completion_metadata: dict[int, tuple[str | None, object | None]] = {}
        self.pending: dict[int, tuple[Future, SupportsRunT, bool]] = {}
        self.scheduler_shutdown = scheduler_shutdown
        self.checksum_pool = ThreadPoolExecutor(
            max_workers=min(checksum_workers, _CACHE_WRITER_CHECKSUM_THREADS),
            thread_name_prefix="haddock-cache-checksum",
        )
        self._thread = threading.Thread(
            target=self._run,
            name="haddock-cache-writer",
            daemon=True,
        )

    def start(self) -> None:
        """Start regular cache-record cycles."""
        self._thread.start()

    def join_after_scheduler_shutdown(self) -> None:
        """Join after the scheduler has marked itself shut down."""
        self._thread.join()

    def _run(self) -> None:
        while not self.scheduler_shutdown.is_set():
            self.cycle()
            self.scheduler_shutdown.wait(_CACHE_WRITER_POLL_SECONDS)

    def _drain_completion_queue(self) -> None:
        while True:
            try:
                completion = self.completion_queue.get_nowait()
            except queue.Empty:
                return
            if isinstance(completion, tuple):
                identifier, checksum, mapping = completion
                self.completion_metadata[identifier] = (checksum, mapping)
            else:
                identifier = completion
            self.completed.add(identifier)

    def _collect_pending(self, wait: bool = False) -> None:
        """Append finished records; this writer thread is the sole appender."""
        for identifier, (future, task, split) in tuple(self.pending.items()):
            if not wait and not future.done():
                continue
            try:
                payload = future.result()
                if split:
                    task.append_cache_success_record(payload)
                self.outstanding.remove(identifier)
            except Exception as error:
                # Keep failed checksum jobs outstanding. The next regular cycle
                # retries them, while the final cycle classifies any remainder.
                log.warning("Could not write CNS cache record: %s", error)
            finally:
                self.pending.pop(identifier)

    def _submit_success_record(self, identifier: int, task: SupportsRunT) -> None:
        """Submit independent normalization/checksum work without appending."""
        prepare = getattr(task, "prepare_cache_success_record", None)
        if prepare is None:
            # Compatibility for small scheduler test tasks. Production CNS jobs
            # use the split prepare/append protocol above.
            future = self.checksum_pool.submit(task.write_cache_success_record)
            self.pending[identifier] = (future, task, False)
        else:
            checksum, mapping = self.completion_metadata.get(identifier, (None, None))
            future = self.checksum_pool.submit(prepare, checksum, mapping)
            self.pending[identifier] = (future, task, True)

    def cycle(self, final: bool = False) -> None:
        """Record visible complete jobs, or classify leftovers on final export."""
        self._drain_completion_queue()
        self._collect_pending()
        for identifier in tuple(self.outstanding):
            if identifier in self.pending:
                continue
            task = self.tasks[identifier]
            try:
                # A completion signal is required even on the final cycle.
                # `cache_outputs_present` only asks whether the declared files
                # exist, and a job killed between publishing its PDB and its
                # PSF -- or one whose artifact a previous run left behind --
                # satisfies that while never having produced a result.
                # Recording it would enter a non-result into the store as a
                # valid answer, which is the one failure a cache must never
                # have.  Anything still outstanding is classified as FAILED
                # just below, which is the honest verdict for a job whose
                # completion nobody witnessed.
                if task.cache_outputs_present() and identifier in self.completed:
                    self._submit_success_record(identifier, task)
            except Exception as error:
                # Cache bookkeeping must never bring down a docking run.  Keep
                # the job outstanding so the final cycle can record failure.
                log.warning("Could not write CNS cache record: %s", error)

        if not final:
            return

        self._collect_pending(wait=True)

        for identifier in tuple(self.outstanding):
            task = self.tasks[identifier]
            try:
                task.write_cache_failure_record()
                self.outstanding.remove(identifier)
            except Exception as error:
                log.warning("Could not write FAILED CNS cache record: %s", error)

        self.checksum_pool.shutdown(wait=True)

    def cancel_pending(self) -> None:
        """Stop cache calculations during scheduler interruption without appending."""
        for future, _task, _split in self.pending.values():
            future.cancel()
        self.checksum_pool.shutdown(wait=False, cancel_futures=True)

    def flush_success_records(self) -> None:
        """Wait for submitted success checksums without classifying failures."""
        self._collect_pending(wait=True)


def split_tasks(lst: Sequence[AnyT], n: int) -> Generator[Sequence[AnyT], None, None]:
    """Split tasks into N-sized chunks."""
    # An empty list would make `n` zero below, and `range(0, 0, 0)` raises
    # ValueError. Nothing to yield in that case.
    if not lst:
        return
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

    def __init__(
        self,
        tasks: Sequence[SupportsRunT],
        results: Queue,
        cache_completion_queue: Queue | None = None,
        task_indices: Sequence[int] | None = None,
    ) -> None:
        super(Worker, self).__init__()
        self.tasks = tasks
        self.result_queue = results
        self.cache_completion_queue = cache_completion_queue
        self.task_indices = list(task_indices or range(len(tasks)))
        log.debug(f"Worker ready with {len(self.tasks)} tasks")

    def run(self) -> None:
        """Execute tasks."""
        results = []
        for task in self.tasks:
            r = None
            try:
                r = task.run()
            except CachedCNSFailure:
                # A trusted FAILED cache record is an expected terminal state:
                # CNS is deliberately skipped and the missing output is
                # classified by the cache writer without a worker warning.
                pass
            except ConfigurationError as error:
                # Configuration policy errors must reach the scheduler.  Letting
                # the child die here would leave the parent blocked on its queue.
                self.result_queue.put(error)
                return
            except Exception as e:
                # A failing task yields ``None`` and the worker carries on, so
                # one bad job does not discard its siblings' results.  The
                # worker cannot vanish here: it still reaches its completion
                # marker below.
                log.warning(f"Exception in task execution: {e}")
            finally:
                if (
                    self.cache_completion_queue is not None
                    and hasattr(task, "cache_writer_id")
                ):
                    completion = getattr(task, "cache_writer_completion", None)
                    self.cache_completion_queue.put(
                        completion() if completion is not None else task.cache_writer_id
                    )

            results.append(r)

        # Put results into the queue
        self.result_queue.put((self.task_indices, results))

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
        cache_context=None,
        cache_debug: bool = False,
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
        self.cache_completion_queue: Queue = Queue()
        self.results: list = []
        self.is_shutdown = threading.Event()

        for task in tasks:
            if hasattr(task, "cache_context"):
                task.cache_context = cache_context
                task.cache_debug = cache_debug

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

        cache_aware_tasks = [
            task
            for task in sorted_task_list
            if (
                getattr(task, "cache_context", None) is not None
                and hasattr(task, "cache_outputs_present")
            )
        ]
        for identifier, task in enumerate(cache_aware_tasks):
            task.cache_writer_id = identifier
        manifests, invariant_checksums = precompute_invariant_checksums_for_jobs(
            cache_aware_tasks
        )
        for task in cache_aware_tasks:
            if hasattr(task, "input_file"):
                task.cache_invariant_checksums = invariant_checksums
        for step_dir, dependencies in manifests.items():
            write_cns_dependencies(step_dir, tuple(dependencies))
        self.cache_writer = (
            CacheRecordWriter(
                cache_aware_tasks,
                self.cache_completion_queue,
                self.is_shutdown,
                self.num_processes,
            )
            if cache_aware_tasks
            else None
        )

        self.task_batches = self._prioritized_task_batches(
            sorted_task_list, cache_context
        )
        self.worker_list = self._make_workers(self.task_batches[0])

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
            if self.cache_writer is not None:
                self.cache_writer.start()
            results_by_index: dict = {}
            for task_batch in self.task_batches:
                self.worker_list = self._make_workers(task_batch)
                for worker in self.worker_list:
                    worker.start()
                results_by_index.update(self._collect_worker_results())
                for worker in self.worker_list:
                    worker.join()

            if self.cache_writer is not None:
                # Mark shutdown and join before the final regular cycle so
                # only one thread can remove a job from the outstanding set.
                self.is_shutdown.set()
                self.cache_writer.join_after_scheduler_shutdown()
                # A worker can finish just after the writer's previous poll.
                # Consume those completion signals now; unresolved jobs remain
                # for export's final cycle.
                self.cache_writer.cycle()
                self.cache_writer.flush_success_records()

            # Workers finish independently, and cache candidates are run in a
            # batch of their own, so neither queue arrival order nor batch
            # order is submission order.  Several callers pair these values
            # with their submitted jobs positionally.
            self.results = [
                results_by_index[index] for index in sorted(results_by_index)
            ]

            log.info(f"{self.num_tasks} tasks finished")

        except KeyboardInterrupt as err:
            # Q: why have a keyboard interrupt here?
            # A: To have a controlled break if the user Ctrl+c during CNS run
            self.terminate()
            # this raises sends the error to libs.libworkflow.Step
            # if Scheduler is used independently the error will propagate to
            # whichever has to catch it
            raise err

    def _prioritized_task_batches(self, tasks, cache_context):
        """Run jobs with a source-cache output artifact before likely misses."""
        indexed = list(enumerate(tasks))
        if not getattr(cache_context, "source_indexes", ()):
            return [indexed]

        candidates = []
        misses = []
        for entry in indexed:
            candidate = getattr(entry[1], "has_cached_output_file", lambda: False)
            (candidates if candidate() else misses).append(entry)
        if not candidates:
            return [indexed]
        return [candidates, misses] if misses else [candidates]

    def _make_workers(self, indexed_tasks):
        """Create balanced workers for one scheduling batch.

        Each worker is handed the submission index of every task it runs, so
        results can be restored to submission order however the batches were
        prioritised and however the workers interleave.
        """
        return [
            Worker(
                [task for _, task in jobs],
                self.queue,
                self.cache_completion_queue,
                [index for index, _ in jobs],
            )
            for jobs in split_tasks(indexed_tasks, self.num_processes)
        ]

    def _collect_worker_results(self) -> dict:
        """Collect one batch's results, keyed by submission index."""
        results_by_index: dict = {}
        completed_workers = 0
        while completed_workers < len(self.worker_list):
            result = self.queue.get()
            if isinstance(result, Exception):
                raise result
            if isinstance(result, str) and result.endswith("_done"):
                completed_workers += 1
            else:
                task_indices, results = result
                results_by_index.update(zip(task_indices, results))
        return results_by_index


    def finalize_cache_records(self) -> None:
        """Run the final cache cycle immediately before module IO export."""
        if self.cache_writer is not None:
            self.cache_writer.cycle(final=True)

    def terminate(self) -> None:
        """Terminate tasks in a controlled way."""
        self.is_shutdown.set()
        for worker in self.worker_list:
            worker.terminate()

        if self.cache_writer is not None and self.cache_writer._thread.is_alive():
            self.cache_writer.join_after_scheduler_shutdown()
        if self.cache_writer is not None:
            self.cache_writer.cancel_pending()

        log.info("The workers terminated in a controlled way")
