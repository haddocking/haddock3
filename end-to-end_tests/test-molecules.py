"""End-to-end test for the refine-molecules example which is used to test supported molecules."""

import gzip
import shutil
import tempfile
import warnings

from pathlib import Path
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from haddock.clis.cli import main as cli_main


EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "refine-molecules"
)


def test_refine_molecules(monkeypatch):
    """Test supported molecules (refine-molecules example).

    Uses refine-molecules.cfg to run the refinement pipeline
    (topoaa, mdref)
    """
    warnings.simplefilter("ignore", PDBConstructionWarning)

    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        config_fname = "refine-molecules.cfg"
        cfg = Path(tmpdir, config_fname)
        shutil.copy(Path(EXAMPLE_DIR, config_fname), cfg)

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-test")

        # Verify all workflow steps produced output directories
        assert Path(run_dir, "0_topoaa").exists(), f"0_topoaa not created"
        assert Path(run_dir, "1_mdref").exists(), f"1_mdref not created"
        assert Path(run_dir, "1_mdref/mdref_1.pdb.gz").exists(), f"1_mdref/mdref_1.pdb.gz not created"