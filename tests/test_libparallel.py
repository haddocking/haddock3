import uuid
from multiprocessing import Queue
from pathlib import Path

import pytest

from haddock.libs.libparallel import (
    GenericTask,
    Scheduler,
    Worker,
    get_index_list,
    split_tasks,
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
