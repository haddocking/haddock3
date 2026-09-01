"""Test libhpc."""

import os
import pytest
import pytest_mock  # noqa : F401

from pathlib import Path
from subprocess import CompletedProcess

from haddock.libs.libhpc import (
    HPCScheduler,
    HPCWorker,
    extract_slurm_status,
    JOB_STATUS_DIC,
    to_torque_time,
)

from haddock.libs.libsubprocess import CNSJob
from haddock.core.exceptions import CNSRunningError


def test_to_torque_time():
    """Test minutes to HH:MM:SS cast."""
    assert to_torque_time(10) == "00:10:00"
    assert to_torque_time(60) == "01:00:00"
    assert to_torque_time(70) == "01:10:00"
    assert to_torque_time(1510) == "25:10:00"
    assert to_torque_time(6070) == "101:10:00"


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


def test_slurm_status(slurm_scontrol_terminated_jobid):
    status = extract_slurm_status(slurm_scontrol_terminated_jobid)
    assert status == "RUNNING"
    assert JOB_STATUS_DIC[status] == "running"


@pytest.fixture
def slurm_scontrol_wrongjobid():
    return "slurm_load_jobs error: Invalid job id specified"


def test_slurm_nojobid(slurm_scontrol_wrongjobid):
    status = extract_slurm_status(slurm_scontrol_wrongjobid)
    assert status == "FAILED"


@pytest.fixture
def hpcworker(mocker):
    """Instanciate a HPCWorker object."""
    mocker.patch(
        "haddock.libs.libsubprocess.CNSJob.cns_exec",
        return_value=None,
    )
    return HPCWorker(
        tasks=[
            CNSJob(
                Path("rigidbody.inp"),
                Path("rigidbody.out"),
                envvars={
                    "MODDIR": ".",
                    "TOPPAR": "topology_params",
                    "MODULE": "rigidbody",
                },
                cns_exec=None,
            ),
        ],
        num=1,
        job_id=123456789,
        workfload_manager="slurm",
        queue=10,
    )


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
    os.remove(hpcworker.job_fname)
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


def test_hpcworker_normalizes_task_outputs(tmp_path, monkeypatch):
    """Test that HPC workers normalize CNS output artifacts after execution."""
    work_dir = tmp_path / "1_rigidbody"
    work_dir.mkdir()
    (work_dir / "rigidbody_1.pdb").write_text("REMARK DATE: volatile\nATOM\n")
    cns_exec = tmp_path / "cns"
    cns_exec.write_text("#!/bin/sh\n")
    cns_exec.chmod(0o755)

    monkeypatch.chdir(work_dir)
    job = CNSJob(
        Path("rigidbody_1.inp"),
        Path("rigidbody_1.out"),
        envvars={
            "MODDIR": str(work_dir),
            "TOPPAR": "topology_params",
            "MODULE": "rigidbody",
        },
        cns_exec=cns_exec,
        output_files=[Path("rigidbody_1.pdb")],
    )
    worker = HPCWorker(
        tasks=[job],
        num=1,
        job_id=123456789,
        workfload_manager="slurm",
    )
    monkeypatch.chdir(tmp_path)

    worker.normalize_outputs()

    assert (work_dir / "rigidbody_1.pdb").read_text() == "ATOM\n"


def test_hpcworker_publishes_tasks_independently(hpcworker, mocker, tmp_path):
    """One faulty concatenated task must not discard another task's output."""
    failed_task = hpcworker.tasks[0]
    successful_task = mocker.MagicMock()
    failed_publish = mocker.patch.object(
        failed_task,
        "publish_outputs",
        side_effect=CNSRunningError("missing output"),
    )
    hpcworker.tasks.append(successful_task)
    partial_input = tmp_path / "1_1.partial.inp"
    partial_input.write_text("input\n")
    hpcworker.partial_inputs = [partial_input]

    hpcworker.normalize_outputs()

    failed_publish.assert_called_once_with(check_output_log=True)
    successful_task.publish_outputs.assert_called_once_with(check_output_log=True)
    assert not partial_input.exists()


def test_hpcworker_job_file_continues_after_each_concatenated_task(
    tmp_path,
    monkeypatch,
):
    """A hard failure in one shell command must not skip later CNS tasks."""
    work_dir = tmp_path / "1_rigidbody"
    work_dir.mkdir()
    executable = tmp_path / "cns"
    executable.write_text("#!/bin/sh\n")
    executable.chmod(0o755)
    monkeypatch.chdir(work_dir)
    jobs = []
    for index in (1, 2):
        input_file = Path(f"rigidbody_{index}.inp")
        input_file.write_text(f'eval ($output_pdb_filename="rigidbody_{index}.pdb")\n')
        jobs.append(
            CNSJob(
                input_file,
                Path(f"rigidbody_{index}.out"),
                envvars={
                    "MODDIR": str(work_dir),
                    "TOPPAR": "topology_params",
                    "MODULE": "rigidbody",
                },
                cns_exec=executable,
                output_files=[Path(f"rigidbody_{index}.pdb")],
            )
        )
    worker = HPCWorker(tasks=jobs, num=1, workfload_manager="slurm")

    worker.prepare_job_file()

    job_file = worker.job_fname.read_text()
    assert "|| exit" not in job_file
    assert f"{executable} < 1_1.partial.inp > rigidbody_1.out" in job_file
    assert f"{executable} < 1_2.partial.inp > rigidbody_2.out" in job_file


def test_hpcscheduler_normalizes_failed_workers(hpcworker, mocker):
    """Terminal worker status must not suppress its successful task outputs."""
    scheduler = HPCScheduler(hpcworker.tasks)
    worker = scheduler.worker_list[0]
    mocker.patch.object(worker, "run")

    def mark_failed():
        worker.job_status = "failed"
        return worker.job_status

    mocker.patch.object(worker, "update_status", side_effect=mark_failed)
    normalize = mocker.patch.object(worker, "normalize_outputs")

    scheduler.run()

    normalize.assert_called_once_with()
