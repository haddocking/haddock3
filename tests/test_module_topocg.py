"""Specific tests for topocg."""

import random
import shutil
import tempfile
import warnings
from math import isnan
from pathlib import Path
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import pytest

import haddock.modules.topology.topocg as topocg_mod
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs import libpdb
from haddock.libs.libaa2cg import add_dummy
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


def test_dummy_bead_placement_uses_the_supplied_seed():
    beads = [("SC1", [1.0, 2.0, 3.0])]

    first = add_dummy(beads, rng=random.Random(494))
    second = add_dummy(beads, rng=random.Random(494))
    different = add_dummy(beads, rng=random.Random(495))

    assert first == second
    assert first != different


def test_run_uses_presplit_models_without_splitting(monkeypatch, protein):
    """Models from the previous (topoaa) step are already split.

    ``_run`` must consume them directly and must not re-split the ensemble
    (the redundant ``libpdb.split_ensemble`` call that was removed).
    """
    with tempfile.TemporaryDirectory() as tempdir:
        root = Path(tempdir)
        prev_dir = root / "1_topoaa"
        cwd = root / "2_topocg"
        prev_dir.mkdir()
        cwd.mkdir()

        # Two individual (already-split) models with their aa topologies,
        # mimicking the output of the previous topoaa step.
        prev_output = {}
        expected_inputs = []
        for idx in range(2):
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

        # Guard against regressions: splitting must never happen here.
        def _no_split(*args, **kwargs):
            raise AssertionError("topocg should not call split_ensemble")

        monkeypatch.setattr(libpdb, "split_ensemble", _no_split)

        # Capture the model paths handed to topology generation instead of
        # running CNS, and avoid the DSSP/martinize machinery.
        captured_inputs = []

        def _fake_generate_topology(input_pdb, *args, **kwargs):
            captured_inputs.append(Path(input_pdb))
            return Path(f"{Path(input_pdb).stem}.{Format.CNS_INPUT}")

        monkeypatch.setattr(topocg_mod, "generate_topology", _fake_generate_topology)

        # Neutralise the CNS engine.
        class _FakeEngine:
            def __init__(self, jobs):
                self.jobs = jobs

            def run(self):
                return None

        monkeypatch.setattr(topocg_mod, "get_engine", lambda *a, **k: _FakeEngine)

        # Skip the disk-dependent export step, just capture its input.
        exported = {}

        def _fake_export(self, faulty_tolerance=0.0):
            exported["models"] = self.output_models

        monkeypatch.setattr(Topocg, "export_io_models", _fake_export)

        module._run()

    # The presplit model paths were used verbatim (no split, no renaming).
    assert captured_inputs == expected_inputs
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
