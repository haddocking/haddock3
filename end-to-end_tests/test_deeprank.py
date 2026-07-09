"""End-to-end test for the deeprank scoring example."""

import shutil
import tempfile
import warnings

from pathlib import Path

from Bio.PDB.PDBExceptions import PDBConstructionWarning

from haddock.clis.cli import main as cli_main
from tests.conftest import has_deeprank

EXAMPLE_DIR = Path(__file__).resolve().parents[1] / "examples" / "scoring"


@has_deeprank
def test_deeprank_scoring(monkeypatch):
    """Test deeprank scoring example.

    Uses deeprank-test.cfg to run the full pipeline (topoaa, deeprank,
    caprieval) and checks that the scoring tsv file is produced.
    """
    warnings.simplefilter("ignore", PDBConstructionWarning)

    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        config_fname = "deeprank-test.cfg"
        cfg = Path(tmpdir, config_fname)
        shutil.copy(Path(EXAMPLE_DIR, config_fname), cfg)

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-deeprank-test")

        assert Path(run_dir, "0_topoaa").exists(), "0_topoaa not created"
        assert Path(run_dir, "1_deeprank").exists(), "1_deeprank not created"
        assert Path(run_dir, "1_deeprank", "deeprank.tsv").exists(), (
            "1_deeprank/deeprank.tsv not created"
        )
        assert Path(run_dir, "2_caprieval").exists(), "2_caprieval not created"
