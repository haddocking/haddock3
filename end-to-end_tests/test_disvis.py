"""End-to-end test for the disvis module."""

import shutil
from pathlib import Path

import pytest

from haddock.clis.cli import main as cli_main
from haddock.libs.libplots import read_capri_table
from haddock.modules.sampling.disvis.disvis import disvis_is_available


if not disvis_is_available():
    pytest.skip("disvis>=3.3.0 dependency not found", allow_module_level=True)

EXAMPLE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "docking-protein-protein"
)

CFG_CONTENT = """
run_dir = "run1-disvis-test"
ncores = 2
molecules = ["data/e2aP_1F3G.pdb", "data/hpr_1.pdb"]

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

[disvis]
ambig_fname = "data/e2a-hpr_air.tbl"
angle = 20
sampling = 2

[flexref]
ambig_fname = "data/e2a-hpr_air.tbl"

[caprieval]
reference_fname = "data/e2a-hpr_1GGR.pdb"
"""


def test_topoaa_disvis_flexref(monkeypatch, tmp_path):
    """The disvis models are read and refined by flexref."""
    data = Path(tmp_path, "data")
    data.mkdir()
    for fname in ("e2aP_1F3G.pdb", "e2a-hpr_air.tbl", "e2a-hpr_1GGR.pdb"):
        shutil.copy(Path(EXAMPLE_DIR, "data", fname), data)
    # The first model of the HPR ensemble
    ensemble = Path(EXAMPLE_DIR, "data", "hpr_ensemble.pdb").read_text()
    Path(data, "hpr_1.pdb").write_text(ensemble.split("ENDMDL")[0].split("\n", 1)[1])
    Path(tmp_path, "disvis-test.cfg").write_text(CFG_CONTENT)

    monkeypatch.chdir(tmp_path)
    cli_main(Path("disvis-test.cfg"))

    capri = read_capri_table(Path("run1-disvis-test", "3_caprieval", "capri_ss.tsv"))
    assert len(capri) == 2
    assert all(model.startswith("../2_flexref/") for model in capri["model"])
