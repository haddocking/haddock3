"""End-to-end test for the protein-ligand docking example."""

import shutil
import tempfile
from pathlib import Path

from haddock.clis.cli import main as cli_main

EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "docking-protein-ligand"
)
INTEGRATION_DIR = (
    Path(__file__).resolve().parents[1] / "integration_tests" / "golden_data"
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
        assert Path(run_dir, "04_caprieval").exists()
        assert Path(run_dir, "05_flexref").exists()
        assert Path(run_dir, "06_caprieval").exists()
        assert Path(run_dir, "07_ilrmsdmatrix").exists()
        assert Path(run_dir, "08_clustrmsd").exists()
        assert Path(run_dir, "09_seletopclusts").exists()
        assert Path(run_dir, "10_caprieval").exists()
        assert Path(run_dir, "11_caprieval").exists()


def test_protein_ligand_autotoppar_ens_workflow(monkeypatch):
    workflow = """
run_dir = "autotoppar-ens"
molecules = ["neuraminidase-2BAT.pdb", "oseltamivir_ens.pdb"]
debug = true
[topoaa]
autohis = true
delenph = false
autotoppar = true
[rigidbody]
tolerance = 0
ambig_fname = "ambig-active-rigidbody.tbl"
w_vdw = 1.0
sampling = 3
[seletop]
select = 1
[flexref]
tolerance = 0
ambig_fname = "ambig-passive.tbl"
mdsteps_rigid = 0
mdsteps_cool1 = 0
"""

    receptor = EXAMPLE_DIR / "data" / "neuraminidase-2BAT.pdb"
    ligand = INTEGRATION_DIR / "oseltamivir_ens.pdb"
    restraint_rigidbody = EXAMPLE_DIR / "data" / "ambig-active-rigidbody.tbl"
    restraint = EXAMPLE_DIR / "data" / "ambig-passive.tbl"
    input_files = [receptor, ligand, restraint_rigidbody, restraint]

    with tempfile.TemporaryDirectory(delete=False) as tmpdir:
        for src in input_files:
            dst = Path(tmpdir, src.name)
            shutil.copy(src, dst)

        with open(f"{tmpdir}/run.toml", "w") as fh:
            fh.write(workflow)

        monkeypatch.chdir(tmpdir)

        print(f"tmpdir: {tmpdir}")
        cli_main("run.toml")
