"""Specific tests for topocg."""

import random
import shutil
import tempfile
import math
from math import isnan
from pathlib import Path
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import pytest

import haddock.modules.topology.topocg as topocg_mod
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs import libpdb
from haddock.libs.libaa2cg import DEFAULT_SEED, add_dummy, martinize
from haddock.libs.libontology import Format, PDBFile
from haddock.modules.topology.topocg import DEFAULT_CONFIG as topocg_params
from haddock.modules.topology.topocg import HaddockModule as Topocg
from haddock.modules.topology.topocg import generate_topology

from . import golden_data


DEFAULT_DICT = read_from_yaml_config(topocg_params)


@pytest.mark.parametrize(
    "param",
    ["hisd_1", "hise_1"],
)
def test_variable_defaults_are_nan_in_mol1(param):
    """Test some variable defaults are as expected."""
    assert isnan(DEFAULT_DICT["mol1"][param])


def test_default_seed_matches_the_topocg_default():
    """``libaa2cg`` hardcodes the fallback seed; keep it equal to the schema.

    Callers without an ``iniseed`` of their own (``caprieval``, ``caprifilter``)
    fall back to it, and should coarse-grain a reference the same way ``topocg``
    coarse-grains the models.
    """
    assert DEFAULT_SEED == DEFAULT_DICT["iniseed"]


def test_there_is_only_one_mol():
    """Test there is only one mol parameter in topocg."""
    r = set(p for p in DEFAULT_DICT if p.startswith("mol") and p[3].isdigit())
    assert len(r) == 1


@pytest.fixture(name="ensemble_header_w_md5")
def fixture_ensemble_header_w_md5():
    """???"""
    return Path(golden_data, "ens_header.pdb")


@pytest.fixture(name="protein")
def fixture_protein():
    """???"""
    return Path(golden_data, "protein.pdb")


@pytest.fixture(name="topocg")
def fixture_topocg(monkeypatch):
    """topocg module fixture"""
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        yield Topocg(order=1, path=Path("."), initial_params=topocg_params)


def test_generate_topology(topocg, protein):
    """Test generate_topology function."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        force_field = topocg.params["cgffversion"]
        observed_inp_out = generate_topology(
            input_pdb=protein,
            output_path=topocg.path,
            recipe_str=topocg.recipe_str,
            defaults=topocg.params,
            mol_params=topocg.params.pop("mol1"),
            default_params_path=None,
            force_field=force_field,
            seed=topocg.params["iniseed"],
        )

    assert observed_inp_out == Path(protein.name).with_suffix(f".{Format.CNS_INPUT}")

    # Confirm the expected output file exists
    expected_filename = topocg.path.resolve() / f"{protein.stem}_cg.pdb"
    assert expected_filename.exists(), (
        f"Expected CG PDB file not found: {expected_filename}"
    )

    # Check that backmapping tbl file has been created
    tbl_backmapping_fpath = topocg.path.resolve() / f"{protein.stem}_cg_to_aa.tbl"
    assert tbl_backmapping_fpath.exists(), (
        f"Expected backmapping retraints file not found: {tbl_backmapping_fpath}"
    )

    # Inspect the content to confirm it's a CG structure
    with open(expected_filename, "r") as f:
        contents = f.read()

    # Check for CG-specific atom/residue names (e.g., "BB" for backbone bead)
    assert "BB" in contents, "CG markers not found in output PDB."


class _FakeEngine:
    """Stand-in for the CNS engine: accepts the jobs, runs nothing."""

    def __init__(self, jobs):
        self.jobs = jobs

    def run(self):
        return None


def _prepared_topocg(monkeypatch, tempdir, protein, n_models):
    """Wire up a topocg module over ``n_models`` presplit inputs, without CNS.

    Returns ``(module, expected_inputs, captured)``, where ``captured`` collects
    the keyword arguments handed to ``generate_topology`` for each model.
    """
    root = Path(tempdir)
    prev_dir = root / "1_topoaa"
    cwd = root / "2_topocg"
    prev_dir.mkdir()
    cwd.mkdir()

    # Individual (already-split) models with their aa topologies, mimicking
    # the output of the previous topoaa step.
    prev_output = {}
    expected_inputs = []
    for idx in range(n_models):
        pdb_name = f"mol_{idx + 1}.pdb"
        shutil.copy(protein, prev_dir / pdb_name)
        (prev_dir / f"mol_{idx + 1}.{Format.TOPOLOGY}").write_text("")
        pdbfile = PDBFile(file_name=pdb_name, path=str(prev_dir))
        prev_output[idx] = pdbfile
        expected_inputs.append(pdbfile.rel_path)

    monkeypatch.chdir(cwd)
    module = Topocg(order=1, path=Path("."), initial_params=topocg_params)
    module.envvars = {}

    class _FakePreviousIO:
        output = [prev_output]

    module.previous_io = _FakePreviousIO()

    # Capture what topology generation is asked for instead of running CNS,
    # which also keeps the DSSP/martinize machinery out of these tests.
    captured = []

    def _fake_generate_topology(input_pdb, *args, **kwargs):
        captured.append({"input_pdb": Path(input_pdb), **kwargs})
        return Path(f"{Path(input_pdb).stem}.{Format.CNS_INPUT}")

    monkeypatch.setattr(topocg_mod, "generate_topology", _fake_generate_topology)
    monkeypatch.setattr(topocg_mod, "get_engine", lambda *a, **k: _FakeEngine)
    return module, expected_inputs, captured


def test_run_uses_presplit_models_without_splitting(monkeypatch, protein):
    """Models from the previous (topoaa) step are already split.

    ``_run`` must consume them directly and must not re-split the ensemble
    (the redundant ``libpdb.split_ensemble`` call that was removed).
    """
    with tempfile.TemporaryDirectory() as tempdir:
        module, expected_inputs, captured = _prepared_topocg(
            monkeypatch, tempdir, protein, n_models=2
        )

        # Guard against regressions: splitting must never happen here.
        def _no_split(*args, **kwargs):
            raise AssertionError("topocg should not call split_ensemble")

        monkeypatch.setattr(libpdb, "split_ensemble", _no_split)

        # Skip the disk-dependent export step, just capture its input.
        exported = {}

        def _fake_export(self, faulty_tolerance=0.0):
            exported["models"] = self.output_models

        monkeypatch.setattr(Topocg, "export_io_models", _fake_export)

        module._run()

    # The presplit model paths were used verbatim (no split, no renaming).
    assert [c["input_pdb"] for c in captured] == expected_inputs
    # One CNS job per input model was created.
    assert len(module.output_models) == 1
    assert len(module.output_models[0]) == 2
    assert exported["models"] is module.output_models


def test_get_md5(topocg, ensemble_header_w_md5, protein):
    """Test get_md5 method."""
    observed_md5_dic = topocg.get_md5(ensemble_header_w_md5)
    expected_md5_dic = {
        1: "71098743056e0b95fbfafff690703761",
        2: "f7ab0b7c751adf44de0f25f53cfee50b",
        3: "41e028d8d28b8d97148dc5e548672142",
        4: "761cb5da81d83971c2aae2f0b857ca1e",
        5: "6c438f941cec7c6dc092c8e48e5b1c10",
    }

    assert observed_md5_dic == expected_md5_dic

    observed_md5_dic = topocg.get_md5(protein)
    assert observed_md5_dic == {}


def _cg_pdb_text(pdb, output_dir, **kwargs):
    """Coarse-grain ``pdb`` into ``output_dir`` and return the CG PDB text."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        cg_pdb = martinize(pdb, str(output_dir), True, **kwargs)
    return Path(cg_pdb).read_text()


def test_martinize_is_reproducible(protein):
    """Same seed must give the same CG model, run after run.

    The ``SCD*`` dummy beads are placed along a random vector; an unseeded
    generator made every ``topocg`` run produce different CG coordinates.
    """
    with tempfile.TemporaryDirectory() as tempdir:
        first = _cg_pdb_text(protein, Path(tempdir, "a"))
        second = _cg_pdb_text(protein, Path(tempdir, "b"))

    assert first == second


def test_martinize_seed_changes_dummy_beads(protein):
    """A different seed must reorient the dummy beads."""
    with tempfile.TemporaryDirectory() as tempdir:
        first = _cg_pdb_text(protein, Path(tempdir, "a"), seed=1)
        other = _cg_pdb_text(protein, Path(tempdir, "b"), seed=2)

    assert first != other
    # Only the dummy beads may move: everything else is deterministic mapping.
    differing = {
        line[12:16].strip()
        for line, other_line in zip(first.splitlines(), other.splitlines())
        if line.startswith("ATOM") and line != other_line
    }
    assert differing == {"SCD1", "SCD2"}


def test_martinize_does_not_touch_global_random(protein):
    """``martinize`` must use its own generator, not the global one.

    Otherwise the result would depend on whatever else consumed the global
    ``random`` stream earlier in the same interpreter.
    """
    random.seed(42)
    expected = [random.random(), random.random()]

    random.seed(42)
    observed = [random.random()]
    with tempfile.TemporaryDirectory() as tempdir:
        _cg_pdb_text(protein, Path(tempdir))
    # If martinize drew from the global stream, the next value would have moved.
    observed.append(random.random())

    assert observed == expected


def test_martinize_is_order_independent(protein, ensemble_header_w_md5):
    """Converting another structure first must not change the result."""
    with tempfile.TemporaryDirectory() as tempdir:
        alone = _cg_pdb_text(protein, Path(tempdir, "a"))
        _cg_pdb_text(ensemble_header_w_md5, Path(tempdir, "b"))
        after = _cg_pdb_text(protein, Path(tempdir, "c"))

    assert alone == after


def test_run_passes_iniseed_to_generate_topology(monkeypatch, protein):
    """``_run`` must forward the ``iniseed`` parameter to the CG conversion."""
    with tempfile.TemporaryDirectory() as tempdir:
        module, _, captured = _prepared_topocg(
            monkeypatch, tempdir, protein, n_models=1
        )
        module.params["iniseed"] = 4242
        monkeypatch.setattr(Topocg, "export_io_models", lambda self, **k: None)

        module._run()

    assert captured[0]["seed"] == 4242


def test_add_dummy_matches_force_field_distances():
    """Dummy beads must be placed at the ``SCd`` bond lengths, in angstrom.

    ``cns/toppar/protein-CG-Martini-2-2.param`` defines ``BOND * SCd = 1.100``
    and ``BOND SCd SCd = 2.800``, so a pair at +/-1.4 A straddling the parent
    is 2.8 A apart and a lone bead sits 1.1 A away.
    """
    rng = random.Random(0)
    parent = (0.0, 0.0, 0.0)

    pair = add_dummy([("SC1", parent)], rng, dist=1.4, n=2)
    assert math.dist(parent, pair["SCD1"]) == pytest.approx(1.4)
    assert math.dist(parent, pair["SCD2"]) == pytest.approx(1.4)
    # The pair must straddle the parent, not sit on the same side of it.
    assert math.dist(pair["SCD1"], pair["SCD2"]) == pytest.approx(2.8)

    single = add_dummy([("SC1", parent)], rng, dist=1.1, n=1)
    assert list(single) == ["SCD1"]
    assert math.dist(parent, single["SCD1"]) == pytest.approx(1.1)


def test_add_dummy_directions_are_isotropic():
    """Bead orientations must be uniform on the sphere, not cube-biased.

    Normalising a vector sampled from a cube pulls directions towards the
    cube's corners and away from the axes. For a uniform distribution the
    component along any axis is uniform on [-1, 1], so a fraction 0.1 of the
    directions satisfies ``|z| > 0.9``; cube sampling gives roughly 0.06.
    """
    rng = random.Random(1)
    n_samples = 20000
    parent = (0.0, 0.0, 0.0)

    near_axis = 0
    for _ in range(n_samples):
        bead = add_dummy([("SC1", parent)], rng, dist=1.0, n=1)
        if abs(bead["SCD1"][2]) > 0.9:
            near_axis += 1

    # ~5 sigma for n_samples=20000; cube sampling (~0.06) fails this.
    assert near_axis / n_samples == pytest.approx(0.1, abs=0.01)
