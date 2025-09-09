import time
from pathlib import Path

import pytest

from haddock.libs.libgrid import GridJob, JobStatus, ping_dirac

from . import has_grid


@pytest.fixture
def grid_job():
    instructions = """#!/bin/bash
echo 'Hello World'
hostname
touch output.txt

zip output.zip output.txt
exit 0
"""
    yield GridJob(instructions)


@has_grid
def test_ping_dirac():
    assert ping_dirac()


@has_grid
def test_gridjob_run(grid_job):

    assert grid_job.status == JobStatus.UNKNOWN
    assert grid_job.id is None

    grid_job.run()

    assert grid_job.status != JobStatus.UNKNOWN
    assert grid_job.id is not None


@has_grid
def test_gridjob_status(grid_job):

    grid_job.run()

    for _ in range(60):
        grid_job.update_status()
        print(grid_job)
        time.sleep(0.5)
        if grid_job.status == JobStatus.DONE:
            break

    assert grid_job.status == JobStatus.DONE


@has_grid
def test_gridjob_retrieve_output(grid_job):
    grid_job.id = 211987118
    output_list = grid_job.retrieve_output()

    assert len(output_list) == 1
    assert Path(output_list[0]).name == "output.txt"
