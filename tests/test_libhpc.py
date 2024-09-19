"""Test libhpc."""
import os
import pytest
import pytest_mock  # noqa : F401

from pathlib import Path
from subprocess import CompletedProcess

from haddock.libs.libhpc import (
    HPCWorker,
    extract_slurm_status,
    JOB_STATUS_DIC,
    to_torque_time,
    )

from haddock.libs.libsubprocess import CNSJob


def test_to_torque_time():
    """Test minutes to HH:MM:SS cast."""
    assert to_torque_time(10) == '00:10:00'
    assert to_torque_time(60) == '01:00:00'
    assert to_torque_time(70) == '01:10:00'
    assert to_torque_time(1510) == '25:10:00'
    assert to_torque_time(6070) == '101:10:00'


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
    assert status == 'RUNNING'
    assert JOB_STATUS_DIC[status] == 'running'


@pytest.fixture
def slurm_scontrol_wrongjobid():
    return 'slurm_load_jobs error: Invalid job id specified'


def test_slurm_nojobid(slurm_scontrol_wrongjobid):
    status = extract_slurm_status(slurm_scontrol_wrongjobid)
    assert status == 'FAILED'


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
                Path('rigidbody.inp'),
                Path('rigidbody.out'),
                envvars={
                    'MODDIR': '.',
                    'TOPPAR': 'topology_params',
                    'MODULE': 'rigidbody',
                    },
                cns_exec=None,
                ),
            ],
        num=1,
        job_id=123456789,
        workfload_manager='slurm',
        queue=10,
        )


def test_hpcworker_run(hpcworker, mocker):
    """Test the `run` function of a HPCWorker object."""
    mocker.patch(
        "subprocess.run",
        return_value=CompletedProcess(
            args=['sbatch', str(hpcworker.job_fname)],
            returncode=0,
            stdout=b'Submitted batch job 42914957',
            stderr=b'',
            )
        )
    hpcworker.run()
    assert os.path.exists(hpcworker.job_fname)
    os.remove(hpcworker.job_fname)
    assert hpcworker.job_id == 42914957
    assert hpcworker.job_status == 'submitted'


def test_hpcworker_update_status(
    hpcworker,
    slurm_scontrol_terminated_jobid,
    mocker,
        ):
    """Test `update_status` function."""
    mocker.patch(
        "subprocess.run",
        return_value=CompletedProcess(
            args=['scontrol', 'show', 'jobid', '-dd', '42909924'],
            returncode=0,
            stdout=bytes(slurm_scontrol_terminated_jobid, 'utf-8'),
            stderr=b'',
            )
        )
    status = hpcworker.update_status()
    assert status == hpcworker.job_status
    assert status == 'running'
