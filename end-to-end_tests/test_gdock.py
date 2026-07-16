"""End-to-end test for the gdock implementation."""

import gzip
import shutil
import tempfile
import warnings

from pathlib import Path
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from haddock.clis.cli import main as cli_main
from haddock.libs.libplots import read_capri_table


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

[flexref]
ambig_fname = "data/e2a-hpr_air.tbl"

[caprieval]
reference_fname = "data/e2a-hpr_1GGR.pdb"
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


def test_topoaa_gdock_flexref(monkeypatch):
    """Test the topoaa -> gdock -> flexref interaction.

    Verifies that the complex PDB produced by gdock is correctly passed on
    to flexref (i.e. flexref can read its coordinates/topology and produce a
    refined model).
    """
    warnings.simplefilter("ignore", PDBConstructionWarning)

    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(Path(EXAMPLE_DIR, "data"), Path(tmpdir, "data"))
        _write_single_model(
            Path(tmpdir, "data", "hpr_ensemble.pdb"),
            Path(tmpdir, "data", "hpr_ensemble_1.pdb"),
        )

        cfg = Path(tmpdir, "docking-protein-protein-gdock-flexref-test.cfg")
        cfg.write_text(CFG_CONTENT)

        monkeypatch.chdir(tmpdir)
        cli_main(cfg)

        run_dir = Path("run1-gdock-test")

        assert Path(run_dir, "0_topoaa").exists(), "0_topoaa not created"
        assert Path(run_dir, "1_gdock").exists(), "1_gdock not created"
        assert Path(run_dir, "2_flexref").exists(), "2_flexref not created"
        assert Path(run_dir, "3_caprieval").exists(), "3_caprieval not created"

        flexref_pdb = Path(run_dir, "2_flexref", "flexref_1.pdb.gz")
        assert flexref_pdb.exists(), "flexref_1.pdb.gz not created"

        with gzip.open(flexref_pdb, "rt") as f:
            chains = {line[21] for line in f if line.startswith(("ATOM", "HETATM"))}
        assert chains == {"A", "B"}, f"Expected chains A and B, got {chains}"

        # Check structural integrity via caprieval metrics
        capri_ss = Path(run_dir, "3_caprieval", "capri_ss.tsv")
        assert capri_ss.exists(), "capri_ss.tsv not created by caprieval"

        df = read_capri_table(capri_ss)
        assert len(df) > 0, "capri_ss.tsv is empty"

        row = df.iloc[0]

        # FIXME: This is not really checking anything
        assert row["irmsd"] >= 0.0, "irmsd must be non-negative"
        assert 0.0 <= row["fnat"] <= 1.0, "fnat must be in [0, 1]"
        assert 0.0 <= row["dockq"] <= 1.0, "dockq must be in [0, 1]"
