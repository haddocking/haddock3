from haddock.clis import cli_score
from tests import golden_data
from pathlib import Path
import tempfile
import pytest_mock
import os
import shutil
import io
from contextlib import redirect_stdout


def test_cli_score_main(mocker):
    """Test the main function of the cli_score module."""
    # setup
    pdb_f = Path(golden_data, "protprot_complex_1.pdb")
    tmpdir = "tmpdir"
    os.mkdir(tmpdir)
    os.mkdir(Path(tmpdir, "1_emscoring"))
    shutil.copy(pdb_f, Path(tmpdir, "1_emscoring", "emscoring_1.pdb"))
    # mocking
    mocker.patch("haddock.libs.libworkflow.WorkflowManager.run", return_value=None)
    mocker.patch("shutil.rmtree", return_value=None)
    mocker.patch("pathlib.Path.mkdir", return_value=None)
    # parsing
    f = io.StringIO()
    with redirect_stdout(f):
        cli_score.main(pdb_f, tmpdir, full=True, keep_all=True)
    out = f.getvalue().split(os.linesep)
    # clean up
    os.unlink(Path(tmpdir, "1_emscoring", "emscoring_1.pdb"))
    os.rmdir(Path(tmpdir, "1_emscoring"))
    os.rmdir(tmpdir)
