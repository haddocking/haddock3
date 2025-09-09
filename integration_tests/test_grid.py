from haddock.libs.libgrid import GridJob
import pytest
from haddock.libs.libsubprocess import CNSJob
import time

from haddock.libs.libgrid import JobStatus

from . import CNS_EXEC


@pytest.fixture
def cnsjob():
    """A simple CNSJob that writes a message to a file."""
    expected_output_fname = "output.pdb"
    cns_inp = f"""
set display={expected_output_fname} end
display hello_from_the_grid
close {expected_output_fname} end
stop"""

    yield CNSJob(
        input_file=cns_inp,
        output_file=expected_output_fname,
        cns_exec=CNS_EXEC,
    )


def test_submit_cns_job_to_grid(cnsjob):
    instructions = """#!/bin/bash
unzip payload.zip
./cns < input.inp 
"""
    j = GridJob(cnsjob, instructions)
    j.submit()

    # HACK: WAIT FOREVER
    while j.status != JobStatus.DONE:
        j.update_status()
        print(j)
        time.sleep(2)
        if j.status == JobStatus.FAILED:
            pytest.fail("Job failed on the grid")

    j.retrieve_output()
