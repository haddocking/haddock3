import pytest

from haddock.libs.libgrid import (
    GridJob,
    CompositeGridJob,
    JobStatus,
    Tag,
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


def test_tag_enum():
    assert Tag.PROBING.value == "Probing"
    assert Tag.DEFAULT.value == "Default"


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
