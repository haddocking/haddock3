"""End-to-end test for the explicit membrane protein-protein docking example."""

import gzip
import shutil
import tempfile
from pathlib import Path

from haddock.clis.cli import main as cli_main

EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "docking-membrane-protein"
)


def test_protein_ligand_autotoppar_workflow(monkeypatch):
    """Test explicit membrane protein-protein docking example.

    Uses docking-membrane-protein-one-side-test.cfg to run the full pipeline
    (topoaa, rigidbody, caprieval, seletop, flexref, caprieval, emref
    caprieval, clustfcc, seletopclusts, caprieval`.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        cfg = Path(tmpdir, "docking-membrane-protein-one-side-test.cfg")
        shutil.copy(
            Path(EXAMPLE_DIR, "docking-membrane-protein-one-side-test.cfg"), cfg
        )

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-test-mem-one-side")

        # Verify all workflow steps produced output directories
        assert Path(run_dir, "00_topoaa").exists()
        assert Path(run_dir, "01_rigidbody").exists()
        assert Path(run_dir, "02_caprieval").exists()
        assert Path(run_dir, "03_seletop").exists()
        assert Path(run_dir, "04_flexref").exists()
        assert Path(run_dir, "05_caprieval").exists()
        assert Path(run_dir, "06_emref").exists()
        assert Path(run_dir, "07_caprieval").exists()
        assert Path(run_dir, "08_clustfcc").exists()
        assert Path(run_dir, "09_seletopclusts").exists()
        assert Path(run_dir, "10_caprieval").exists()

        emref_file = Path(run_dir, "06_emref", "emref_1.pdb.gz")
        with gzip.open(emref_file, "rt") as f:
            content = f.read()
            assert "DPP" in content, "DPP string not found in emref_1.pdb.gz"
