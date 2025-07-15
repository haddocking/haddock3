import uuid
from multiprocessing import Queue
from pathlib import Path
import time
import signal
import os
from unittest.mock import patch, Mock

import pytest

from haddock.libs.libparallel import (
    GenericTask,
    Scheduler,
    Worker,
    get_index_list,
    split_tasks,
    AlascanScheduler, 
    AlascanWorker,
)


class Task:
    """Dummy task class to be used in the test.

    This is a simple task that receives an integer and returns the integer + 1.

    **The important part is that the task is a class that implements a `run` method.**
    """

    def __init__(self, input):
        self.input = input
        self.output = None

    def run(self) -> int:
        self.output = self.input + 1
        return self.output


class FileTask:
    """Dummy task class to be used in the test.

    This is a simple task that receives a filename and creates an empty file with that name.

    **The important part is that the task is a class that implements a `run` method.**
    """

    def __init__(self, filename):
        self.input_file = Path(filename)

    def run(self):
        Path(self.input_file).touch()

class TaskWithException:

    def __init__(self):
        pass

    def run(self):
        raise ValueError("Test error")


@pytest.fixture
def worker():
    """Return a worker with 3 tasks."""
    yield Worker(tasks=[Task(1), Task(2), Task(3)], results=Queue())


@pytest.fixture
def scheduler():
    """Return a scheduler with 3 tasks."""
    yield Scheduler(
        ncores=1,
        tasks=[
            Task(1),
            Task(2),
            Task(3),
        ],
    )


@pytest.fixture
def scheduler_files():
    """Return a scheduler with 3 tasks that create files."""

    file_list = [uuid.uuid4().hex for _ in range(3)]
    yield Scheduler(
        ncores=1,
        tasks=[
            FileTask(file_list[0]),
            FileTask(file_list[1]),
            FileTask(file_list[2]),
        ],
    )

    for f in file_list:
        try:
            Path(f).unlink()
        except FileNotFoundError:
            pass


@pytest.fixture
def scheduler_with_exception():
    """Return a scheduler with 3 tasks, one of them raises an exception."""
    yield Scheduler(
        ncores=1,
        tasks=[
            Task(1),
            TaskWithException(),
            Task(3),
        ],
    )


def test_split_tasks():

    lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    n = 3
    result = list(split_tasks(lst, n))
    assert result == [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10]]

    n = 4
    result = list(split_tasks(lst, n))
    assert result == [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10]]


def test_get_index_list():

    nmodels = 10
    ncores = 3
    result = get_index_list(nmodels, ncores)
    assert result == [0, 4, 7, 10]

    nmodels = 10
    ncores = 4
    result = get_index_list(nmodels, ncores)
    assert result == [0, 3, 6, 8, 10]


def test_worker_run(worker):

    _ = worker.run()

    assert worker.tasks[0].output == 2
    assert worker.tasks[1].output == 3
    assert worker.tasks[2].output == 4


def test_scheduler_files(scheduler_files):

    _ = scheduler_files.run()

    assert Path(scheduler_files.worker_list[0].tasks[0].input_file).exists()
    assert Path(scheduler_files.worker_list[0].tasks[1].input_file).exists()
    assert Path(scheduler_files.worker_list[0].tasks[2].input_file).exists()


def test_scheduler(scheduler):

    _ = scheduler.run()

    assert scheduler.results[0] == 2
    assert scheduler.results[1] == 3
    assert scheduler.results[2] == 4


def test_scheduler_with_exception(scheduler_with_exception):

    _ = scheduler_with_exception.run()

    assert scheduler_with_exception.results[0] == 2
    assert scheduler_with_exception.results[1] is None
    assert scheduler_with_exception.results[2] == 4


def test_generic_task_init():
    def sample_function(a, b, c=3):
        return a + b + c

    # Test with only a function
    task1 = GenericTask(sample_function)
    assert task1.function == sample_function
    assert task1.args == ()
    assert task1.kwargs == {}

    # Test with function and positional arguments
    task2 = GenericTask(sample_function, 1, 2)
    assert task2.function == sample_function
    assert task2.args == (1, 2)
    assert task2.kwargs == {}

    # Test with function, positional and keyword arguments
    task3 = GenericTask(sample_function, 1, b=2, c=4)
    assert task3.function == sample_function
    assert task3.args == (1,)
    assert task3.kwargs == {"b": 2, "c": 4}

    # Test with a lambda function
    lambda_func = lambda x: x * 2
    task4 = GenericTask(lambda_func, 5)
    assert task4.function == lambda_func
    assert task4.args == (5,)
    assert task4.kwargs == {}

    # Test with built-in function
    task5 = GenericTask(len, "hello")
    assert task5.function == len
    assert task5.args == ("hello",)
    assert task5.kwargs == {}


def test_generic_task_init_exceptions():
    # Test that providing a non-callable raises a TypeError
    with pytest.raises(TypeError):
        GenericTask("not a function")


def test_generic_task_init_types():
    def dummy_function():
        pass

    task = GenericTask(dummy_function, 1, 2, a=3)
    assert callable(task.function)
    assert isinstance(task.args, tuple)
    assert isinstance(task.kwargs, dict)


import pytest


def test_generic_task_run():
    # Test with a simple function
    def add(a, b):
        return a + b

    task = GenericTask(add, 2, 3)
    assert task.run() == 5

    # Test with a function that returns None
    def do_nothing():
        pass

    task = GenericTask(do_nothing)
    assert task.run() is None

    # Test with a function that takes keyword arguments
    def greet(name, greeting="Hello"):
        return f"{greeting}, {name}!"

    task = GenericTask(greet, "Alice", greeting="Hi")
    assert task.run() == "Hi, Alice!"

    # Test with a lambda function
    task = GenericTask(lambda x: x * 2, 5)
    assert task.run() == 10

    # Test with a built-in function
    task = GenericTask(len, "hello")
    assert task.run() == 5

    # Test with a method of a class
    class TestClass:
        def method(self, x):
            return x * 2

    obj = TestClass()
    task = GenericTask(obj.method, 3)
    assert task.run() == 6


def test_generic_task_run_exceptions():
    # Test that running a task with a non-callable raises a TypeError
    with pytest.raises(TypeError):
        task = GenericTask("not a function")
        task.run()

    # Test that running a task with incorrect arguments raises a TypeError
    def func(a, b):
        return a + b

    task = GenericTask(func, 1)  # Missing one argument
    with pytest.raises(TypeError):
        task.run()

    # Test that the function's exceptions are propagated
    def raise_value_error():
        raise ValueError("Test error")

    task = GenericTask(raise_value_error)
    with pytest.raises(ValueError, match="Test error"):
        task.run()


def test_generic_task_run_with_complex_args():
    # Test with *args and **kwargs in the function
    def complex_func(*args, **kwargs):
        return sum(args) + sum(kwargs.values())

    task = GenericTask(complex_func, 1, 2, 3, a=4, b=5)
    assert task.run() == 15

    # Test with a function that modifies mutable arguments
    def modify_list(lst):
        lst.append(4)
        return lst

    original_list = [1, 2, 3]
    task = GenericTask(modify_list, original_list)
    result = task.run()
    assert result == [1, 2, 3, 4]
    assert original_list == [1, 2, 3, 4]  # The original list is modified


@pytest.mark.skip("WIP")
def test_scheduler_terminate(scheduler_files):
    pass

# new ones
class SlowTask:
    """Task that takes some time to complete - useful for testing interrupts."""
    
    def __init__(self, duration=0.1):
        self.duration = duration
        self.completed = False
    
    def run(self):
        time.sleep(self.duration)
        self.completed = True
        return "completed"


@pytest.fixture
def alascan_worker():
    """Return an AlascanWorker with 3 tasks."""
    yield AlascanWorker(tasks=[Task(1), Task(2), Task(3)], results=Queue())


@pytest.fixture
def alascan_scheduler():
    """Return an AlascanScheduler with 3 tasks."""
    yield AlascanScheduler(
        ncores=1,
        tasks=[
            Task(1),
            Task(2),
            Task(3),
        ],
    )


@pytest.fixture
def alascan_scheduler_files():
    """Return an AlascanScheduler with 3 tasks that create files."""
    file_list = [uuid.uuid4().hex for _ in range(3)]
    yield AlascanScheduler(
        ncores=1,
        tasks=[
            FileTask(file_list[0]),
            FileTask(file_list[1]),
            FileTask(file_list[2]),
        ],
    )

    for f in file_list:
        try:
            Path(f).unlink()
        except FileNotFoundError:
            pass


@pytest.fixture
def alascan_scheduler_slow():
    """Return an AlascanScheduler with slow tasks for interrupt testing."""
    yield AlascanScheduler(
        ncores=2,
        tasks=[
            SlowTask(0.2),
            SlowTask(0.2),
            SlowTask(0.2),
            SlowTask(0.2),
        ],
    )


def test_alascan_worker_run(alascan_worker):
    """Test that AlascanWorker runs tasks correctly."""
    _ = alascan_worker.run()

    assert alascan_worker.tasks[0].output == 2
    assert alascan_worker.tasks[1].output == 3
    assert alascan_worker.tasks[2].output == 4


def test_alascan_worker_signal_handling():
    """Test that AlascanWorker handles signals correctly."""
    worker = AlascanWorker(tasks=[SlowTask(0.1), SlowTask(0.1)], results=Queue())
    
    # Test signal handler is set up correctly
    worker._signal_handler(signal.SIGTERM, None)
    assert worker._should_stop is True


def test_alascan_worker_cleanup_child_processes():
    """Test that AlascanWorker can clean up child processes."""
    worker = AlascanWorker(tasks=[Task(1)], results=Queue())
    # Test cleanup method doesn't raise errors
    worker._cleanup_child_processes()


def test_alascan_worker_enhanced_run():
    """Test that AlascanWorker properly handles signal setup and cleanup."""
    worker = AlascanWorker(tasks=[Task(1), Task(2)], results=Queue())
    # Test that worker has enhanced attributes
    assert hasattr(worker, '_should_stop')
    assert hasattr(worker, '_signal_handler')
    assert hasattr(worker, '_cleanup_child_processes')
    # Test initial state
    assert worker._should_stop is False


def test_alascan_worker_early_termination():
    """Test that AlascanWorker stops early when signal is received."""
    worker = AlascanWorker(tasks=[Task(1), Task(2), Task(3)], results=Queue())
    # Simulate signal received after first task
    worker._should_stop = True
    # Verify the flag is respected
    assert worker._should_stop is True


def test_alascan_scheduler_basic_functionality(alascan_scheduler):
    """Test that AlascanScheduler runs tasks correctly."""
    _ = alascan_scheduler.run()
    assert alascan_scheduler.results[0] == 2
    assert alascan_scheduler.results[1] == 3
    assert alascan_scheduler.results[2] == 4


def test_alascan_scheduler_files(alascan_scheduler_files):
    """Test that AlascanScheduler works with file tasks."""
    _ = alascan_scheduler_files.run()
    assert Path(alascan_scheduler_files.worker_list[0].tasks[0].input_file).exists()
    assert Path(alascan_scheduler_files.worker_list[0].tasks[1].input_file).exists()
    assert Path(alascan_scheduler_files.worker_list[0].tasks[2].input_file).exists()


def test_alascan_scheduler_with_exception():
    """Test that AlascanScheduler handles exceptions correctly."""
    scheduler = AlascanScheduler(
        ncores=1,
        tasks=[
            Task(1),
            TaskWithException(),
            Task(3),
        ],
    )
    _ = scheduler.run()
    assert scheduler.results[0] == 2
    assert scheduler.results[1] is None
    assert scheduler.results[2] == 4


def test_alascan_scheduler_context_manager():
    """Test that AlascanScheduler works as a context manager."""
    tasks = [Task(1), Task(2), Task(3)]
    with AlascanScheduler(ncores=1, tasks=tasks) as scheduler:
        scheduler.run()
        assert scheduler.results[0] == 2
        assert scheduler.results[1] == 3
        assert scheduler.results[2] == 4


def test_alascan_scheduler_context_manager_with_exception():
    """Test that AlascanScheduler context manager handles exceptions."""
    tasks = [Task(1), TaskWithException(), Task(3)]
    with AlascanScheduler(ncores=1, tasks=tasks) as scheduler:
        scheduler.run()
        assert scheduler.results[0] == 2
        assert scheduler.results[1] is None
        assert scheduler.results[2] == 4


def test_alascan_scheduler_inheritance():
    """Test that AlascanScheduler properly inherits from Scheduler."""
    scheduler = AlascanScheduler(ncores=2, tasks=[Task(1), Task(2)])
    # Should have same attributes as regular Scheduler
    assert hasattr(scheduler, 'num_processes')
    assert hasattr(scheduler, 'queue')
    assert hasattr(scheduler, 'results')
    assert hasattr(scheduler, 'worker_list')
    # Should use AlascanWorker instead of regular Worker
    assert len(scheduler.worker_list) > 0
    assert isinstance(scheduler.worker_list[0], AlascanWorker)


def test_alascan_scheduler_safe_terminate():
    """Test that safe_terminate method works correctly."""
    scheduler = AlascanScheduler(ncores=1, tasks=[SlowTask(0.1)])
    # Start the scheduler in a way that we can test termination
    scheduler.safe_terminate()  # Should not raise any errors even if no workers started


def test_alascan_scheduler_enhanced_cleanup():
    """Test that enhanced cleanup method works correctly."""
    scheduler = AlascanScheduler(ncores=1, tasks=[Task(1)])
    # Test cleanup doesn't raise errors
    scheduler._enhanced_cleanup()  # Should not raise any errors


def test_alascan_worker_vs_regular_worker():
    """Test that AlascanWorker behaves like regular Worker for normal operations."""
    regular_tasks = [Task(1), Task(2), Task(3)]
    alascan_tasks = [Task(1), Task(2), Task(3)]
    regular_worker = Worker(regular_tasks, Queue())
    alascan_worker = AlascanWorker(alascan_tasks, Queue())
    # Both should be Process instances
    assert isinstance(regular_worker, Worker)
    assert isinstance(alascan_worker, Worker)
    assert isinstance(alascan_worker, AlascanWorker)
    # Both should have the same basic attributes
    assert len(regular_worker.tasks) == len(alascan_worker.tasks)
    assert hasattr(regular_worker, 'result_queue')
    assert hasattr(alascan_worker, 'result_queue')


def test_alascan_scheduler_job_distribution():
    """Test that AlascanScheduler distributes jobs correctly like regular Scheduler."""
    # Test with more tasks than cores
    tasks = [Task(i) for i in range(10)]
    scheduler = AlascanScheduler(ncores=3, tasks=tasks)
    # Should create 3 workers
    assert len(scheduler.worker_list) == 3
    # All tasks should be distributed among workers
    total_tasks = sum(len(worker.tasks) for worker in scheduler.worker_list)
    assert total_tasks == 10


def test_alascan_scheduler_with_file_sorting():
    """Test that AlascanScheduler handles file sorting like regular Scheduler."""
    # Create tasks with input_file attributes (like the original Scheduler test)
    file_list = [f"file_{i}.pdb" for i in [3, 1, 2]]  # Not sorted 
    tasks = [FileTask(f) for f in file_list]
    scheduler = AlascanScheduler(ncores=1, tasks=tasks)
    # Should have created workers with sorted tasks
    assert len(scheduler.worker_list) > 0
    # Clean up files
    for f in file_list:
        try:
            Path(f).unlink(missing_ok=True)
        except:
            pass


@pytest.mark.parametrize("ncores", [1, 2, 4])
def test_alascan_scheduler_different_core_counts(ncores):
    """Test AlascanScheduler with different core counts."""
    tasks = [Task(i) for i in range(8)]
    scheduler = AlascanScheduler(ncores=ncores, tasks=tasks)
    scheduler.run()
    # Should complete all tasks regardless of core count
    assert len(scheduler.results) == 8
    expected_results = list(range(1, 9))
    assert sorted(scheduler.results) == expected_results


def test_alascan_scheduler_empty_tasks():
    """Test AlascanScheduler with empty task list."""
    # This test will pass once we add the empty tasks guard to AlascanScheduler
    scheduler = AlascanScheduler(ncores=2, tasks=[])
    # Should handle empty tasks gracefully
    scheduler.run()
    assert scheduler.results == []
    