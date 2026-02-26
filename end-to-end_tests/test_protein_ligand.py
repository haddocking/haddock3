"""End-to-end test for the protein-ligand docking example."""

import shutil
import tempfile
from pathlib import Path

from haddock.clis.cli import main as cli_main

EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "docking-protein-ligand"
)


def test_protein_ligand_autotoppar_workflow(monkeypatch):
    """Test protein-ligand docking example with automated topology generation.

    Uses docking-protein-ligand-autotoppar-test.cfg to run the full pipeline
    (topoaa, rigidbody, caprieval, seletop, flexref, caprieval, rmsdmatrix,
    clustrmsd, seletopclusts, caprieval x2) without any user-provided ligand
    topology or parameter files.  ``autotoppar=True`` instructs topoaa to
    invoke prodrg automatically, and the generated files are propagated to all
    downstream modules via ``_output_params``.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        cfg = Path(tmpdir, "docking-protein-ligand-autotoppar-test.cfg")
        shutil.copy(
            Path(EXAMPLE_DIR, "docking-protein-ligand-autotoppar-test.cfg"), cfg
        )

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-autotoppar-test")

        # Check if the auto-generated prodrg topology files were generated
        autotoppar_param = Path(
            run_dir, "00_topoaa", "oseltamivir_zwitterion_prodrg.param"
        )
        autotoppar_top = Path(run_dir, "00_topoaa", "oseltamivir_zwitterion_prodrg.top")
        assert autotoppar_param.exists(), f"{autotoppar_param} was not generated"
        assert autotoppar_top.exists(), f"{autotoppar_top} was not generated"

        # Verify all workflow steps produced output directories
        assert Path(run_dir, "00_topoaa").exists()
        assert Path(run_dir, "01_rigidbody").exists()
        assert Path(run_dir, "02_caprieval").exists()
        assert Path(run_dir, "03_seletop").exists()
        assert Path(run_dir, "04_flexref").exists()
        assert Path(run_dir, "05_caprieval").exists()
        assert Path(run_dir, "06_rmsdmatrix").exists()
        assert Path(run_dir, "07_clustrmsd").exists()
        assert Path(run_dir, "08_seletopclusts").exists()
        assert Path(run_dir, "09_caprieval").exists()
        assert Path(run_dir, "10_caprieval").exists()
