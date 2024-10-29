import tempfile
from pathlib import Path
import os
import shutil
from haddock.libs.libworkflow import WorkflowManager
from haddock.core.typing import Any
from integration_tests import GOLDEN_DATA


def test_emscoring_workflow(caplog, monkeypatch):
    """Test WorkflowManager."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # copy the pdb file to the temporary directory
        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_1.pdb"),
            Path(tmpdir, "protglyc_complex_1.pdb"),
        )
        monkeypatch.chdir(tmpdir)

        caplog.set_level("INFO")
        ParamDict = {
            'topoaa.1': {
                'autohis': True,
                'molecules': [Path("protglyc_complex_1.pdb")],
                'clean': True,
                'offline': False,
                'mode': "local",
                'ncores': 2,
                },
            'emscoring.1': {
                'per_interface_scoring' : False,
                }
            }

        workflow = WorkflowManager(
            ParamDict,
            start=0,
            other_params=Any,
            )
        workflow.run()
        # assert the directories 0_topoaa and 1_emscoring have been created
        assert os.path.isdir("0_topoaa") is True
        assert os.path.isdir("1_emscoring") is True
        workflow.postprocess()
        # assert an analysis directory has been created
        assert os.path.isdir("analysis") is True
        assert os.path.isdir("analysis/1_emscoring_analysis") is True
        # there should be a report.html inside it
        assert os.path.isfile("analysis/1_emscoring_analysis/report.html") is True
