"""Test libhpc."""

import os
import pytest

from pathlib import Path
from subprocess import CompletedProcess

from haddock.libs.libhpc import (
    HPCWorker,
    extract_slurm_status,
    JOB_STATUS_DIC,
    to_torque_time,
)
import tempfile

from haddock.libs.libsubprocess import CNSJob


@pytest.fixture
def slurm_scontrol_terminated_jobid():
    return """JobId=42909924 JobName=proabc2.job
   UserId=enmr(1095) GroupId=users(100) MCS_label=N/A
   Priority=1 Nice=0 Account=(null) QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   DerivedExitCode=0:0
   RunTime=00:00:01 TimeLimit=04:00:00 TimeMin=N/A
   SubmitTime=2023-10-17T10:03:01 EligibleTime=2023-10-17T10:03:01
   AccrueTime=2023-10-17T10:03:01
   StartTime=2023-10-17T10:03:02 EndTime=2023-10-17T14:03:02 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   LastSchedEval=2023-10-17T10:03:02
   Partition=short AllocNode:Sid=bianca:124111
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=node011
   BatchHost=node011
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,node=1,billing=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
     Nodes=node011 CPU_IDs=4 Mem=0 GRES_IDX=
   MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/trinity/home/enmr/csb_webserver/data/runs/proabc2/nAO_3OnarD4a/proabc2.job
   WorkDir=/trinity/home/enmr/csb_webserver/data/runs/proabc2
   StdErr=/trinity/home/enmr/csb_webserver/data/runs/proabc2/nAO_3OnarD4a/proabc2.err
   StdIn=/dev/null
   StdOut=/trinity/home/enmr/csb_webserver/data/runs/proabc2/nAO_3OnarD4a/proabc2.out
   Power=
"""


@pytest.fixture
def slurm_scontrol_wrongjobid():
    return "slurm_load_jobs error: Invalid job id specified"


@pytest.fixture
def cnsjob(mocker):
    with tempfile.NamedTemporaryFile() as f:
        f.file.write(b"")
        f.file.flush()
        f.file.seek(0)
        os.chmod(f.name, 0o755)
        # mocker.patch(
        #     "haddock.libs.libsubprocess.CNSJob.cns_exec",
        #     return_value=f.name,
        # )
        yield CNSJob(
            Path("rigidbody.inp"),
            Path("rigidbody.out"),
            envvars={
                "MODDIR": ".",
                "TOPPAR": "topology_params",
                "MODULE": "rigidbody",
            },
            cns_exec=f.name,
        )


@pytest.fixture
def hpcworker(cnsjob):
    """Instanciate a HPCWorker object."""
    with tempfile.NamedTemporaryFile() as f:
        _hpcworder = HPCWorker(
            tasks=[cnsjob],
            num=1,
            job_id=123456789,
            workfload_manager="slurm",
            queue=str(10),
        )
        _hpcworder.job_fname = Path(f.name)
        yield _hpcworder


def test_hpcworker_init(cnsjob):
    hpcworker = HPCWorker(
        tasks=[cnsjob],
        num=1,
        job_id=123456789,
        workfload_manager="slurm",
    )

    assert hpcworker.tasks == [cnsjob]
    assert hpcworker.job_num == 1
    assert hpcworker.job_id == 123456789
    assert hpcworker.job_status == "unknown"

    assert hpcworker.moddir == Path(cnsjob.envvars["MODDIR"])
    assert hpcworker.toppar == cnsjob.envvars["TOPPAR"]
    assert hpcworker.cns_folder == cnsjob.envvars["MODULE"]

    module_name = Path(cnsjob.envvars["MODDIR"]).resolve().stem.split("_")[-1]
    assert hpcworker.job_fname == Path(
        Path(cnsjob.envvars["MODDIR"]), f"{module_name}_1.job"
    )

    assert hpcworker.workload_manager == "slurm"

    assert hpcworker.queue is None


def test_hpcworker_prepare_job_file(mocker, hpcworker):

    mocked_create_job_header_funcs = mocker.patch(
        "haddock.libs.libhpc.create_job_header_funcs"
    )
    mocked_create_job_header_funcs.__getitem__.return_value.return_value = (
        f"header{os.linesep}"
    )

    mocker.patch(
        "haddock.libs.libhpc.create_CNS_export_envvars",
        return_value=f"envvars{os.linesep}",
    )

    with tempfile.NamedTemporaryFile() as tmpfile:
        hpcworker.moddir = Path(tmpfile.name)
        hpcworker.prepare_job_file()

        assert hpcworker.job_fname.exists()
        assert hpcworker.job_fname.stat().st_size > 0

        expected_job_contents = f"header{os.linesep}envvars{os.linesep}cd {hpcworker.moddir}{os.linesep}{hpcworker.tasks[0].cns_exec} < {hpcworker.tasks[0].input_file} > {hpcworker.tasks[0].output_file}{os.linesep}"
        with open(hpcworker.job_fname, "r") as f:
            assert f.read() == expected_job_contents


def test_hpcworker_run(hpcworker, mocker):
    """Test the `run` function of a HPCWorker object."""
    mocker.patch(
        "subprocess.run",
        return_value=CompletedProcess(
            args=["sbatch", str(hpcworker.job_fname)],
            returncode=0,
            stdout=b"Submitted batch job 42914957",
            stderr=b"",
        ),
    )
    hpcworker.run()
    assert os.path.exists(hpcworker.job_fname)
    assert hpcworker.job_id == 42914957
    assert hpcworker.job_status == "submitted"


def test_hpcworker_update_status(
    hpcworker,
    slurm_scontrol_terminated_jobid,
    mocker,
):
    """Test `update_status` function."""
    mocker.patch(
        "subprocess.run",
        return_value=CompletedProcess(
            args=["scontrol", "show", "jobid", "-dd", "42909924"],
            returncode=0,
            stdout=bytes(slurm_scontrol_terminated_jobid, "utf-8"),
            stderr=b"",
        ),
    )
    status = hpcworker.update_status()
    assert status == hpcworker.job_status
    assert status == "running"


def test_hpcworker_cancel(mocker, hpcworker):

    mocker.patch.object(hpcworker, "update_status", return_value="running")
    mock_run = mocker.patch(
        "subprocess.run",
        return_value=CompletedProcess(args=[], returncode=0),
    )

    hpcworker.cancel()

    assert hpcworker.job_status == "unknown"

    mock_run.assert_called_once_with(
        ["scancel", str(hpcworker.job_id)],
        capture_output=True,
    )


def test_hpcscheduler_init():
    pass


@pytest.mark.skip(reason="Not implemented yet.")
def test_hpcscheduler_run():
    pass


@pytest.mark.skip(reason="Not implemented yet.")
def test_hpcscheduler_terminate():
    pass


@pytest.mark.skip(reason="Not implemented yet.")
def test_create_slurm_header():
    pass


@pytest.mark.skip(reason="Not implemented yet.")
def test_create_torque_header():
    pass


def test_to_torque_time():
    """Test minutes to HH:MM:SS cast."""
    assert to_torque_time(10) == "00:10:00"
    assert to_torque_time(60) == "01:00:00"
    assert to_torque_time(70) == "01:10:00"
    assert to_torque_time(1510) == "25:10:00"
    assert to_torque_time(6070) == "101:10:00"


def test_slurm_status(slurm_scontrol_terminated_jobid, slurm_scontrol_wrongjobid):
    status = extract_slurm_status(slurm_scontrol_terminated_jobid)
    assert status == "RUNNING"
    assert JOB_STATUS_DIC[status] == "running"

    status = extract_slurm_status(slurm_scontrol_wrongjobid)
    assert status == "FAILED"


@pytest.mark.skip(reason="Not implemented yet.")
def test_create_CNS_export_envvars():
    pass
