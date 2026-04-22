"""End-to-end test for the protein-protein-shape docking example."""

import shutil
import tempfile
import gzip
from pathlib import Path

from haddock.clis.cli import main as cli_main

import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning


EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "docking-protein-protein-shape"
)

def test_protein_protein_shape_CG(monkeypatch):
    """Test protein-protein-shape CG docking example.

    Uses docking-protein-protein-shape-CG-test.cfg to run the full pipeline
    (topoaa, topocg rigidbody, caprieval, seletop, flexref, caprieval, cgtoaa, emref, caprieval,
    clustfcc, seletopclusts, caprieva) 
    """
    warnings.simplefilter('ignore', PDBConstructionWarning)

    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        config_fname = "docking-protein-protein-shape-CG-test.cfg"
        cfg = Path(tmpdir, config_fname)
        shutil.copy(Path(EXAMPLE_DIR, config_fname), cfg)

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-CG-test")

        # Verify all workflow steps produced output directories
        assert Path(run_dir, "00_topoaa").exists(), f"00_topoaa not created"
        assert Path(run_dir, "01_topocg").exists(), f"01_topocg not created"
        assert Path(run_dir, "02_rigidbody").exists(), f"02_rigidbody not created"
        assert Path(run_dir, "03_caprieval").exists(), f"03_caprieval not created"
        assert Path(run_dir, "04_seletop").exists(), f"04_seletop created"
        assert Path(run_dir, "05_flexref").exists(), f"05_flexref not created"
        assert Path(run_dir, "06_caprieval").exists(), f"06_caprieval not created"
        assert Path(run_dir, "07_cgtoaa").exists(), f"07_cgtoaa not created"
        assert Path(run_dir, "08_emref").exists(), f"08_emref not created"
        assert Path(run_dir, "09_caprieval").exists(), f"09_caprieval not created"
        assert Path(run_dir, "10_clustfcc").exists(), f"10_clustfcc not created"
        assert Path(run_dir, "11_seletopclusts").exists(), f"10_seletopclusts not created"
        assert Path(run_dir, "12_caprieval").exists(), f"12_caprieval not created"
        with gzip.open(Path(run_dir, "08_emref/emref_1.pdb.gz"), 'rt') as f:
            assert any('SHA SHA S' in line for line in f), f"Shape atoms not found in emref PDB file"
