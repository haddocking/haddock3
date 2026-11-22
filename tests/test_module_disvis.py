"""Test the disvis module."""

import logging
import shutil
from pathlib import Path

import numpy as np

import pytest

from haddock import log
from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.sampling.disvis import DEFAULT_CONFIG
from haddock.modules.sampling.disvis import HaddockModule as DisvisModule
from haddock.modules.sampling.disvis.disvis import (
    fill_from_clusters,
    pool,
    superposition,
)

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

    def __init__(self, *combinations):
        self.combinations = list(combinations)
        self.output = [m for c in combinations for m in c]

    def retrieve_models(self, crossdock: bool = True):
        """Provide the combinations of molecules."""
        return self.combinations


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
def test_disvis_models(disvis, tmp_path):
    """The models are the receptor and the placed ligand, with hydrogens.

    With an ensemble (two ligand models), `sampling` models are passed on in
    total, from the pooled complexes of both searches, and each search
    leaves the result files of DisVis in its folder.
    """
    receptor = _topoaa_model("e2aP_1F3G_haddock")
    ligand = _topoaa_model("hpr_ensemble_1_haddock")
    other = PDBFile(
        shutil.copy(Path(ligand.path, ligand.file_name), Path(tmp_path, "hpr_2.pdb")),
        topology=ligand.topology,
        path=tmp_path,
    )
    disvis.previous_io = MockPreviousIO((receptor, ligand), (receptor, other))
    disvis._run()

    assert len(disvis.output_models) == 3
    for idx in (1, 2):
        for fname in ("accessible_complexes.out", "violations.out", "args.json"):
            assert Path(f"disvis_search_{idx}", fname).exists()
    natoms = sum(_count_atoms(Path(m.path, m.file_name)) for m in (receptor, ligand))
    for model in disvis.output_models:
        assert _count_atoms(Path(model.file_name)) == natoms
        assert model.topology == [receptor.topology, ligand.topology]


@has_disvis
def test_disvis_two_molecules(disvis):
    """A combination of other than two molecules is an error."""
    receptor = _topoaa_model("e2aP_1F3G_haddock")
    disvis.previous_io = MockPreviousIO((receptor,))
    with pytest.raises(RuntimeError, match="exactly 2 input molecules"):
        disvis._run()


@has_disvis
def test_disvis_log(disvis):
    """The messages of a multiprocessing search go to the haddock3 log, once.

    Logging them on the root logger configured it, which duplicated the
    haddock3 output in the terminal.
    """
    root_handlers = list(logging.root.handlers)
    disvis.params.update(angle=30, ncores=2)
    disvis.previous_io = MockPreviousIO(
        (_topoaa_model("e2aP_1F3G_haddock"), _topoaa_model("hpr_ensemble_1_haddock"))
    )
    records = []
    handler = logging.Handler()
    handler.emit = records.append
    log.addHandler(handler)
    try:
        disvis._run()
    finally:
        log.removeHandler(handler)

    assert logging.root.handlers == root_handlers
    messages = [r.getMessage() for r in records]
    assert sum(m.startswith("Number of processors used") for m in messages) == 1


def test_superposition():
    """A moved receptor model is superposed back on the first one."""
    rng = np.random.default_rng(0)
    reference = rng.normal(0, 10, (20, 3))
    angle = np.radians(40)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    moved = reference @ rotation.T + [5, -3, 8]
    S, s = superposition(moved, reference)
    np.testing.assert_allclose(moved @ S.T + s, reference, atol=1e-10)


@has_disvis
def test_pool_most_consistent_first():
    """The cap keeps the complexes consistent with the most restraints."""
    from disvis.disvis import Complexes

    def part(consistent, key):
        n = len(consistent)
        counts = np.bincount(consistent, minlength=5)
        return Complexes(
            np.arange(n),
            np.zeros((n, 3)),
            np.asarray(consistent),
            np.asarray(key, dtype=np.uint64),
            counts,
        )

    # Search 0 has the complex with 4 restraints, search 1 the two with 3,
    #  one of them with the smallest key of all.
    parts = [part([4, 2, 2], [5, 1, 2]), part([3, 2, 3], [0, 3, 6])]
    kept, search, threshold, nfound = pool(parts, 3, 50, 2)
    assert (threshold, nfound) == (2, 6)
    assert kept.consistent.tolist() == [4, 3]
    assert search.tolist() == [0, 1]
    assert kept.rotation.tolist() == [0, 0]


def test_fill_from_clusters():
    """The centers come first, then other members in rounds over the clusters.

    Within a cluster, the members consistent with the most restraints first.
    """
    clusters = [np.array([0, 1, 2, 3]), np.array([4, 5]), np.array([6, 7])]
    consistent = np.array([5, 3, 4, 4, 5, 4, 5, 4])
    key = np.arange(8)
    assert fill_from_clusters(clusters, consistent, key, 6) == [0, 4, 6, 2, 5, 7]
    assert fill_from_clusters(clusters, consistent, key, 2) == [0, 4]
    assert len(fill_from_clusters(clusters, consistent, key, 20)) == 8
