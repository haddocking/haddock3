"""End-to-end test for the scoring in CG mode example."""

import gzip
import shutil
import tempfile
import warnings

from pathlib import Path
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from haddock.clis.cli import main as cli_main


EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "scoring"
)


def test_protein_protein_shape_CG(monkeypatch):
    """Test coarse-grained scoring example.

    Uses emscoring-cg-test.cfg to run the full pipeline
    (topoaa, topocg, emscoring, cgtoaa, emref) 
    """
    warnings.simplefilter("ignore", PDBConstructionWarning)

    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        config_fname = "emscoring-cg-test.cfg"
        cfg = Path(tmpdir, config_fname)
        shutil.copy(Path(EXAMPLE_DIR, config_fname), cfg)

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-CG-test")

        # Verify all workflow steps produced output directories
        assert Path(run_dir, "0_topoaa").exists(), f"0_topoaa not created"
        assert Path(run_dir, "1_topocg").exists(), f"1_topocg not created"
        assert Path(run_dir, "2_emscoring").exists(), f"2_emscoring not created"
        assert Path(run_dir, "3_cgtoaa").exists(), f"3_cgtoaa not created"
        assert Path(run_dir, "4_emscoring").exists(), f"4_emscoring not created"
        assert Path(run_dir, "4_emscoring/emscoring_1.pdb.gz").exists(), f"04_emscoring/emscoring_1.pdb.gz not created"