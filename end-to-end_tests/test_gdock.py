"""End-to-end test for the topoaa -> gdock interaction."""

import gzip
import shutil
import tempfile
import warnings

from pathlib import Path
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from haddock.clis.cli import main as cli_main


EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "docking-protein-protein"
)


CFG_CONTENT = """
run_dir = "run1-gdock-test"

mode = "local"
ncores = 2

molecules =  [
    "data/e2aP_1F3G.pdb",
    "data/hpr_ensemble_1.pdb"
    ]

[topoaa]
autohis = false
[topoaa.mol1]
nhisd = 0
nhise = 1
hise_1 = 75
[topoaa.mol2]
nhisd = 1
hisd_1 = 76
nhise = 1
hise_1 = 15

[gdock]
ambig_fname = "data/e2a-hpr_air.tbl"
max_generations = 2
sampling = 1
seed = 1
"""


def _write_single_model(src: Path, dst: Path) -> None:
    """Write only the first MODEL of a multi-model PDB file to `dst`."""
    lines = []
    in_model = False
    with open(src) as input_handler:
        for line in input_handler:
            if line.startswith("MODEL"):
                if in_model:
                    break
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                break
            lines.append(line)

    dst.write_text("".join(lines))


def test_topoaa_gdock(monkeypatch):
    """Test the topoaa -> gdock interaction.

    Uses a minimal config (topoaa, gdock) to verify that gdock correctly
    receives the receptor/ligand combinations produced by topoaa.
    """
    warnings.simplefilter("ignore", PDBConstructionWarning)

    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        _write_single_model(
            Path(tmpdir, "data", "hpr_ensemble.pdb"),
            Path(tmpdir, "data", "hpr_ensemble_1.pdb"),
        )

        cfg = Path(tmpdir, "docking-protein-protein-gdock-test.cfg")
        cfg.write_text(CFG_CONTENT)

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-gdock-test")

        assert Path(run_dir, "0_topoaa").exists(), "0_topoaa not created"
        assert Path(run_dir, "1_gdock").exists(), "1_gdock not created"

        gdock_pdb = Path(run_dir, "1_gdock", "gdock_1_1.pdb.gz")
        assert gdock_pdb.exists(), "gdock_1_1.pdb.gz not created"

        with gzip.open(gdock_pdb, "rt") as f:
            chains = {line[21] for line in f if line.startswith(("ATOM", "HETATM"))}
        assert chains == {"A", "B"}, f"Expected chains A and B, got {chains}"
