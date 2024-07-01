import pytest
from haddock.libs.libparallel import split_tasks, get_index_list, Worker, Scheduler
from pathlib import Path
import uuid


class Task:

    def __init__(self, input):
        self.input = input
        self.output = None

    def run(self):
        self.output = self.input + 1


class FileTask:

    def __init__(self, filename):
        self.input_file = Path(filename)

    def run(self):
        Path(self.input_file).touch()


@pytest.fixture
def worker():
    yield Worker(tasks=[Task(1), Task(2), Task(3)])


@pytest.fixture
def scheduler():
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

    # To make sure this was executed, they should have been incremented
    assert worker.tasks[0].output == 2
    assert worker.tasks[1].output == 3
    assert worker.tasks[2].output == 4


def test_scheduler_run(scheduler_files):
    scheduler_files.run()

    assert Path(scheduler_files.worker_list[0].tasks[0].input_file).exists()
    assert Path(scheduler_files.worker_list[0].tasks[1].input_file).exists()
    assert Path(scheduler_files.worker_list[0].tasks[2].input_file).exists()


@pytest.mark.skip("WIP")
def test_scheduler_terminate(scheduler_files):
    pass
