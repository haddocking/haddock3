import tempfile
from pathlib import Path
import os
import shutil

from haddock.clis.cli import main as cli_main
from haddock.clis.cli_analyse import main as cli_analyse
from haddock.clis.cli_re import maincli
from haddock.core.typing import Any
from haddock.libs.libworkflow import WorkflowManager

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


def test_interactive_analysis_on_workflow(monkeypatch):
    """A comprehensive test for the interactive commands and analysis."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # copy the files to the temporary directory
        files = [
            "cyclic-peptide.pdb",
            "prot.pdb",
            "workflow.cfg",
        ]
        for fl in files:
            shutil.copy(
                Path(GOLDEN_DATA, fl),
                Path(tmpdir, fl),
            )

        monkeypatch.chdir(tmpdir)

        
        cli_main(
            Path("workflow.cfg"),
        )
        # read run_dir in workflow.cfg
        with open("workflow.cfg", "r") as f:
            for line in f:
                if "run_dir" in line:
                    run_dir = line.split("=")[1].strip().strip('"')
                    break
        assert os.path.isdir(run_dir) is True
        assert os.path.isdir(Path(run_dir, "0_topoaa")) is True
        assert os.path.isdir(Path(run_dir, "1_rigidbody")) is True
        assert os.path.isdir(Path(run_dir, "2_clustfcc")) is True
        assert os.path.isdir(Path(run_dir, "3_caprieval")) is True
        assert os.path.isdir(Path(run_dir, "analysis")) is True

        # now running interactive re-clustering
        clustfcc_dir = f"{run_dir}/2_clustfcc"
        
        # faking sys.argv in input to haddock3-re
        monkeypatch.setattr("sys.argv",
                            ["haddock3-re", "clustfcc", clustfcc_dir, "-f", "0.7"]
                            )
        maincli()
        assert os.path.isdir(Path(run_dir, "2_clustfcc_interactive")) is True
        assert Path(run_dir, "2_clustfcc_interactive/clustfcc.tsv").exists() is True
        assert Path(run_dir, "2_clustfcc_interactive/clustfcc.txt").exists() is True
        
        # now running interactive re-scoring
        capri_dir = f"{run_dir}/3_caprieval"
        # faking sys.argv in input to haddock3-re
        monkeypatch.setattr("sys.argv",
                            ["haddock3-re", "score", capri_dir, "-a", "1.0"]
                            )
        maincli()
        assert os.path.isdir(Path(run_dir, "3_caprieval_interactive")) is True
        assert Path(run_dir, "3_caprieval_interactive/capri_ss.tsv").exists() is True

        # now analyse the interactive folders
        cli_analyse(
            run_dir,
            [2, 3],
            10,
            format=None,
            scale=None,
            is_cleaned=True,
            inter=True,
            )
        exp_clustfcc_dir = Path(run_dir, "analysis", "2_clustfcc_interactive_analysis")
        exp_caprieval_dir = Path(run_dir, "analysis", "3_caprieval_interactive_analysis")
        assert os.path.isdir(exp_clustfcc_dir) is True
        assert os.path.isdir(exp_caprieval_dir) is True
        assert Path(exp_clustfcc_dir, "report.html").exists() is True
        assert Path(exp_caprieval_dir, "report.html").exists() is True
        
