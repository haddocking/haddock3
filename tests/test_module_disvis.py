"""Test the disvis module."""

import logging
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.sampling.disvis import DEFAULT_CONFIG
from haddock.modules.sampling.disvis import HaddockModule as DisvisModule

from . import golden_data, has_disvis


AMBIG = Path(
    golden_data.parents[1],
    "examples",
    "docking-protein-protein",
    "data",
    "e2a-hpr_air.tbl",
)


class MockPreviousIO:
    """Mocking class to retrieve the topoaa models."""

    def __init__(self, models):
        self.models = models
        self.output = models

    def retrieve_models(self, crossdock: bool = True):
        """Provide the combinations of molecules."""
        return [tuple(self.models)]


def _topoaa_model(name: str) -> PDBFile:
    topology = TopologyFile(Path(golden_data, f"{name}.psf"), path=golden_data)
    return PDBFile(
        Path(golden_data, f"{name}.pdb"), topology=topology, path=golden_data
    )


def _count_atoms(pdb_file: Path) -> int:
    return sum(line.startswith("ATOM") for line in pdb_file.read_text().splitlines())


@pytest.fixture(name="disvis")
def fixture_disvis(monkeypatch, tmp_path):
    """Return a disvis module with a coarse search, in a temporary folder."""
    monkeypatch.chdir(tmp_path)
    module = DisvisModule(order=1, path=Path("."), initial_params=DEFAULT_CONFIG)
    module.params.update(ambig_fname=str(AMBIG), angle=20, sampling=3, ncores=1)
    return module


@has_disvis
def test_disvis_models(disvis):
    """The models are the receptor and the placed ligand, with hydrogens."""
    receptor = _topoaa_model("e2aP_1F3G_haddock")
    ligand = _topoaa_model("hpr_ensemble_1_haddock")
    disvis.previous_io = MockPreviousIO([receptor, ligand])
    disvis._run()

    assert len(disvis.output_models) == 3
    natoms = sum(_count_atoms(Path(m.path, m.file_name)) for m in (receptor, ligand))
    for model in disvis.output_models:
        assert _count_atoms(Path(model.file_name)) == natoms
        assert model.topology == [receptor.topology, ligand.topology]


@has_disvis
def test_disvis_two_molecules(disvis):
    """A combination of other than two molecules is an error."""
    receptor = _topoaa_model("e2aP_1F3G_haddock")
    disvis.previous_io = MockPreviousIO([receptor])
    with pytest.raises(RuntimeError, match="exactly 2 input molecules"):
        disvis._run()


@has_disvis
def test_disvis_log(disvis, caplog):
    """The messages of a multiprocessing search go to the haddock3 log only.

    Logging them on the root logger configured it, which duplicated the
    haddock3 output in the terminal.
    """
    root_handlers = list(logging.root.handlers)
    disvis.params.update(angle=30, ncores=2)
    disvis.previous_io = MockPreviousIO(
        [_topoaa_model("e2aP_1F3G_haddock"), _topoaa_model("hpr_ensemble_1_haddock")]
    )
    disvis._run()

    assert logging.root.handlers == root_handlers
    assert any(
        r.name == "haddock" and r.message.startswith("Number of processors used")
        for r in caplog.records
    )
