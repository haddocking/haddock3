from haddock.libs.libworkflow import WorkflowManager
from haddock.core.typing import Any, ModuleParams
import tempfile
from pathlib import Path
import os

def test_WorkflowManager(caplog):
    """Test WorkflowManager."""
    caplog.set_level("INFO")
    ParamDict = {
        'topoaa.1':
            {'autohis': True,
             'molecules': ['fake.pdb'],
             'clean': True
            }
        }
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        workflow = WorkflowManager(
            ParamDict,
            start=0,
            other_params=Any,)
        workflow.postprocess()
        first_log_line = str(caplog.records[0].message)
        second_log_line = str(caplog.records[1].message)
        assert first_log_line == "Reading instructions step 0_topoaa"
        assert second_log_line == "Running haddock3-analyse on ./, modules [], with top_cluster = 10"
        