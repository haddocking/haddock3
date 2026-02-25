"""End-to-end test for the protein-ligand docking example."""

import shutil
import tempfile
from pathlib import Path

from haddock.clis.cli import main as cli_main

EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"
E2E_DIR = Path(__file__).resolve().parent


def test_protein_ligand_autotoppar_workflow(monkeypatch):
    """Test protein-ligand docking example with automated topology generation.

    Mirrors the docking-protein-ligand-test.cfg pipeline (topoaa, rigidbody,
    seletop, flexref, seletopclusts, caprieval) but without any user-provided
    ligand topology or parameter files.  ``autotoppar=True`` instructs topoaa
    to invoke prodrg automatically, and the generated files are propagated to
    all downstream modules via ``_output_params``.
    """
    example_data = Path(EXAMPLES_DIR, "docking-protein-ligand", "data")

    needed_files = [
        "neuraminidase-2BAT.pdb",
        "oseltamivir_zwitterion.pdb",
        "ambig-active-rigidbody.tbl",
        "ambig-passive.tbl",
        "target.pdb",
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        for fname in needed_files:
            shutil.copy(Path(example_data, fname), Path(tmpdir, fname))

        cfg = Path(tmpdir, "protein_ligand_autotoppar.cfg")
        shutil.copy(Path(E2E_DIR, "protein_ligand_autotoppar.cfg"), cfg)

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run")

        # Auto-generated prodrg topology files (key assertion for autotoppar)
        autotoppar_param = Path(run_dir, "0_topoaa", "oseltamivir_zwitterion_prodrg.param")
        autotoppar_top = Path(run_dir, "0_topoaa", "oseltamivir_zwitterion_prodrg.top")
        assert autotoppar_param.exists(), f"{autotoppar_param} was not generated"
        assert autotoppar_top.exists(), f"{autotoppar_top} was not generated"

        # Verify that all workflow steps produced output directories
        assert Path(run_dir, "0_topoaa").exists()
        assert Path(run_dir, "1_rigidbody").exists()
        assert Path(run_dir, "2_seletop").exists()
        assert Path(run_dir, "3_flexref").exists()
        assert Path(run_dir, "4_caprieval").exists()
