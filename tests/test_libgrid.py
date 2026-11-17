import subprocess
import zipfile
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from haddock.libs.libgrid import (
    GRIDScheduler,
    GridJob,
    CompositeGridJob,
    JobStatus,
    MAX_RETRIES,
    ping_dirac,
    validate_dirac,
)


@pytest.fixture
def dummy_paths(tmp_path):
    toppar = tmp_path / "toppar"
    toppar.mkdir()
    module = tmp_path / "module"
    module.mkdir()
    return toppar, module


def test_gridjob_create_jdl_and_job_script(dummy_paths):
    toppar, module = dummy_paths
    job = GridJob(input="input", toppar_path=toppar, module_path=module)
    job.name = "testjob"
    job.expected_outputs = ["foo.out"]
    job.jdl = job.loc / "job.jdl"
    job.job_script = job.loc / "job.sh"
    job.create_jdl()
    job.create_job_script()
    assert job.jdl.exists()
    assert job.job_script.exists()
    content = job.jdl.read_text()
    assert f'JobName = "{job.name}";' in content
    assert '"foo.out"' in content


def test_gridjob_process_input_f(tmp_path, dummy_paths):
    toppar, module = dummy_paths
    input_str = '$output_test = "result.txt";\n'
    job = GridJob(input=input_str, toppar_path=toppar, module_path=module)
    job.name = "testjob"
    job.loc = tmp_path
    job.input_str = input_str
    job.process_input_f()
    assert "result.txt" in job.expected_outputs
    inp_file = tmp_path / f"{job.name}.inp"
    assert inp_file.exists()


def test_compositegridjob_parse_input_and_process(dummy_paths):
    toppar, module = dummy_paths
    input_list = ["line1\nstop\n", "line2\nstop\n"]
    job = CompositeGridJob(input=input_list, toppar_path=toppar, module_path=module)
    assert len(job.input_str_list) == 2
    job.process_input_f()
    for idx in range(2):
        inp_file = job.loc / f"{idx}_{job.name}.inp"
        assert inp_file.exists()


def test_gridjob_repr(dummy_paths):
    toppar, module = dummy_paths
    job = GridJob(input="input", toppar_path=toppar, module_path=module)
    s = repr(job)
    assert "ID:" in s and "Name:" in s and "Status:" in s


def test_jobstatus_from_string():
    assert JobStatus.from_string("Running") == JobStatus.RUNNING
    assert JobStatus.from_string("unknown") == JobStatus.UNKNOWN
    assert JobStatus.from_string("Stalled") == JobStatus.FAILED


def test_parse_output():
    s = "JobID=123 Status=Running Site=TestSite;"
    d = GridJob.parse_output(s)
    assert d["JobID"] == "123"
    assert d["Status"] == "Running"
    assert d["Site"] == "TestSite"


def test_process_line_and_find_output():
    line = '$output_test = "foo.txt";'
    new_line, found = GridJob._process_line(line)
    assert new_line == line
    output = GridJob._find_output(line)
    assert output == "foo.txt"


def test_clean_timings(dummy_paths):
    toppar, module = dummy_paths
    job = GridJob(input="input", toppar_path=toppar, module_path=module)
    job.timings = {"foo": 123}
    job.clean_timings()
    assert job.timings == {}


def test_clean_removes_tmpdir(dummy_paths):
    toppar, module = dummy_paths
    job = GridJob(input="input", toppar_path=toppar, module_path=module)
    loc = job.loc
    assert loc.exists()
    job.clean()
    assert not loc.exists()


def test_ping_and_validate_dirac(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/" + cmd)
    monkeypatch.setattr(
        "subprocess.run",
        lambda *a, **k: type(
            "R", (), {"returncode": 0, "stderr": b"", "stdout": b""}
        )(),
    )
    assert validate_dirac()
    assert ping_dirac()


def test_process_line(tmp_path, monkeypatch):
    # No match
    line = "no pattern here"
    new_line, found = GridJob._process_line(line)
    assert new_line == line
    assert found is None

    # AT_PATTERN match and file exists
    test_file = tmp_path / "foo.txt"
    test_file.write_text("data")
    line = f"@@{test_file}"
    monkeypatch.chdir(tmp_path)
    new_line, found = GridJob._process_line(line)
    assert new_line == line.replace(str(test_file), "foo.txt")
    assert found == str(test_file)

    # VAR_PATTERN match and file does not exist
    line = '($ test="bar.txt")'
    new_line, found = GridJob._process_line(line)
    assert new_line == line
    assert found is None


def test_find_output():
    line = '$output_test = "foo.txt";'
    output = GridJob._find_output(line)
    assert output == "foo.txt"

    line2 = "$output_test = bar.txt;"
    output2 = GridJob._find_output(line2)
    assert output2 == "bar.txt"

    line3 = "no output here"
    output3 = GridJob._find_output(line3)
    assert output3 is None


@pytest.mark.parametrize(
    "module, n_tasks, expected",
    [("rigidbody", 21, [10, 10, 1]), ("flexref", 1, [1])],
)
def test_scheduler_create_batches(
    tmp_path, monkeypatch, dummy_paths, module, n_tasks, expected
):
    toppar, module_path = dummy_paths
    step = tmp_path / f"1_{module}"
    step.mkdir()
    monkeypatch.chdir(step)
    envvars = {"TOPPAR": toppar, "MODULE": module_path}
    tasks = [
        SimpleNamespace(input_file=f"model{i}\nstop\n", envvars=envvars)
        for i in range(n_tasks)
    ]
    scheduler = GRIDScheduler(tasks, params={"ncores": 1})
    scheduler.create_batches()
    assert [len(job.input_str_list) for job in scheduler.workload] == expected
    for job in scheduler.workload:
        job.clean()


def test_submit_retries_failed_submission(dummy_paths):
    toppar, module = dummy_paths
    job = GridJob(input="input", toppar_path=toppar, module_path=module)
    job.packaged = True
    error = subprocess.CalledProcessError(
        2, "dirac-wms-job-submit", output="Peer closed connection"
    )
    ok = subprocess.CompletedProcess([], 0, stdout="JobID = 123")
    with (
        patch("haddock.libs.libgrid.subprocess.run", side_effect=[error, error, ok]),
        patch("haddock.libs.libgrid.time.sleep") as sleep,
        patch.object(GridJob, "update_status"),
    ):
        job.submit()
    assert job.id == 123
    assert job.retries == 2
    assert [c.args[0] for c in sleep.call_args_list] == [30, 60]
    job.clean()


def test_retrieve_output_skips_missing_cns_log(tmp_path, monkeypatch, dummy_paths):
    toppar, module = dummy_paths
    monkeypatch.chdir(tmp_path)
    job = GridJob(input="input", toppar_path=toppar, module_path=module)
    job.id = 123
    job.expected_outputs = ["model_1.pdb", "0_job.out"]
    sandbox = job.loc / "123"
    sandbox.mkdir()
    (sandbox / "model_1.pdb").write_text("ATOM")
    with patch("haddock.libs.libgrid.subprocess.run"):
        job.retrieve_output()
    assert (tmp_path / "model_1.pdb").read_text() == "ATOM"
    job.clean()


def test_given_up_job_is_cleaned(dummy_paths):
    toppar, module = dummy_paths
    failed, unsubmitted = (
        GridJob(input="input", toppar_path=toppar, module_path=module) for _ in range(2)
    )
    failed.retries = unsubmitted.retries = MAX_RETRIES
    failed.status = JobStatus.FAILED
    unsubmitted.packaged = True
    error = subprocess.CalledProcessError(2, "dirac-wms-job-submit", output="")
    with (
        patch.object(GridJob, "update_status"),
        patch("haddock.libs.libgrid.subprocess.run", side_effect=error),
    ):
        GRIDScheduler.process_job(failed)
        with pytest.raises(subprocess.CalledProcessError):
            unsubmitted.submit()
    assert not failed.loc.exists()
    assert not unsubmitted.loc.exists()


def test_wait_for_completion_checks_all_jobs_at_once(dummy_paths):
    toppar, module = dummy_paths
    jobs = [
        GridJob(input="input", toppar_path=toppar, module_path=module) for _ in range(2)
    ]
    for job_id, job in enumerate(jobs, start=1):
        job.id, job.status = job_id, JobStatus.WAITING
    scheduler = GRIDScheduler([], params={"ncores": 1})
    scheduler.workload = jobs
    status = subprocess.CompletedProcess(
        [],
        0,
        stdout="JobID=1 MinorStatus=Execution Complete; Status=Done; Site=A;\n"
        "JobID=2 MinorStatus=Execution Complete; Status=Done; Site=B;\n",
    )
    with (
        patch("haddock.libs.libgrid.subprocess.run", return_value=status) as run,
        patch("haddock.libs.libgrid.time.sleep") as sleep,
        patch.object(GRIDScheduler, "process_job"),
    ):
        scheduler.wait_for_completion()
    run.assert_called_once()
    assert run.call_args.args[0] == ["dirac-wms-job-status", "1", "2"]
    sleep.assert_called_once_with(30)
    assert [job.site for job in jobs] == ["A", "B"]
    for job in jobs:
        job.clean()


def test_payload_has_no_duplicate_files(dummy_paths):
    toppar, module = dummy_paths
    (toppar / "protein.top").write_text("top")
    job = GridJob(input="input", toppar_path=toppar, module_path=module)
    # The `.inp` refers to the toppar file, so it is also copied next to it
    copied = job.loc / "protein.top"
    copied.write_text("top")
    job.payload_fnames.append(copied)
    job.prepare_payload(cns_script_path=module, toppar_path=toppar)
    with zipfile.ZipFile(job.loc / "payload.zip") as z:
        names = z.namelist()
    assert sorted(names) == ["cns", "protein.top"]
    job.clean()
