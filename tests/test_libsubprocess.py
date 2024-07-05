import pytest
import itertools
import os
from pathlib import Path
import tempfile
import shlex
from unittest.mock import MagicMock
from haddock.libs.libsubprocess import BaseJob, Job, CNSJob


@pytest.fixture
def basejob():
    """Create a BaseJob instance."""
    return BaseJob(
        input_=Path("input"),
        output=Path("output"),
        executable=Path("executable"),
    )


@pytest.fixture
def job():
    """Create a Job instance."""
    return Job(
        input_=Path("input"),
        output=Path("output"),
        executable=Path("executable"),
    )


@pytest.fixture
def cnsjob(mocker):
    """Create a CNSJob instance.

    Here we create a temporary file and set it as the mock CNS executable.
    """
    with tempfile.NamedTemporaryFile() as f:
        f.file.write(b"")
        f.file.flush()
        f.file.seek(0)
        os.chmod(f.name, 0o755)

        mocker.patch("haddock.libs.libsubprocess.global_cns_exec", f.name)

        yield CNSJob(
            input_file=Path("input"),
            output_file=Path("output"),
            cns_exec=Path(f.name),
        )


def test_basejob_run(basejob, mocker):
    """Test the run method of the BaseJob class."""
    basejob.make_cmd = lambda: None
    basejob.cmd = "some_command"

    mock_popen = mocker.patch("subprocess.Popen")
    mock_popen_instance = MagicMock()
    mock_popen.return_value = mock_popen_instance
    mock_popen_instance.communicate.return_value = (b"output", None)

    mock_open = mocker.patch("builtins.open", mocker.mock_open())

    result = basejob.run()

    assert result == b"output"

    basejob.make_cmd()
    assert mock_popen.call_args[0][0] == shlex.split(basejob.cmd)
    mock_open.assert_called_once_with(basejob.output, "w")

    mock_popen.assert_called_once_with(
        shlex.split(basejob.cmd),
        stdout=mock_open(),
        close_fds=True,
    )


def test_basejob_make_cmd(basejob):
    with pytest.raises(NotImplementedError):
        basejob.make_cmd()


def test_jobinputfirst_make_cmd(job):
    job.executable = "executable"
    job.input = "input"
    job.args = [1, 2, 3]

    job.make_cmd()

    assert job.cmd == "executable 123 input"


def test_cnsjob_envvars_setter(cnsjob):
    cnsjob.envvars = {"key": "value"}

    assert cnsjob.envvars == {"key": "value"}
    assert cnsjob._envvars == {"key": "value"}

    with pytest.raises(ValueError):
        cnsjob.envvars = "wrong"


def test_cnsjob_cns_exec_setter(cnsjob, mocker):

    with tempfile.NamedTemporaryFile() as f:
        f.file.write(b"")
        f.file.flush()
        f.file.seek(0)
        os.chmod(f.name, 0o755)

        mocker.patch("haddock.libs.libsubprocess.global_cns_exec", f.name)

        cnsjob.cns_exec = f.name

    with pytest.raises(ValueError):
        cnsjob.cns_exec = "wrong"


def test_cnsjob_run(cnsjob, mocker):

    mock_popen = mocker.patch("subprocess.Popen")
    mock_popen_instance = MagicMock()
    mock_popen.return_value = mock_popen_instance
    mock_popen_instance.communicate.return_value = (b"output", None)

    mocker.patch("builtins.open", mocker.mock_open())
    mocker.patch("haddock.libs.libsubprocess.gzip_files", return_value=None)

    # Try all possible combinations of compress flags
    for comb in itertools.product([True, False], repeat=3):
        compress_inp, compress_out, compress_seed = comb
        result = cnsjob.run(
            compress_inp=compress_inp,
            compress_out=compress_out,
            compress_seed=compress_seed,
        )

        assert result == b"output"
