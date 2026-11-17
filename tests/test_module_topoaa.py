"""Specific tests for topoaa."""

import os
import re
import tempfile
from math import isnan
from pathlib import Path

import pytest

from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.topology.topoaa import DEFAULT_CONFIG as topoaa_params
from haddock.modules.topology.topoaa import HaddockModule as Topoaa
from haddock.modules.topology.topoaa import RECIPE_PATH, generate_topology

from . import golden_data


DEFAULT_DICT = read_from_yaml_config(topoaa_params)


@pytest.mark.parametrize(
    "param",
    ["hisd_1", "hise_1"],
)
def test_variable_defaults_are_nan_in_mol1(param):
    """Test some variable defaults are as expected."""
    assert isnan(DEFAULT_DICT["mol1"][param])


def test_there_is_only_one_mol():
    """Test there is only one mol parameter in topoaa."""
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


@pytest.fixture(name="topoaa")
def fixture_topoaa(monkeypatch):
    """topoaa module fixture"""
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        yield Topoaa(order=1, path=Path("."), initial_params=topoaa_params)


def test_generate_topology(topoaa, protein):
    """Test generate_topology function."""
    observed_inp_out = generate_topology(
        input_pdb=protein,
        recipe_str=topoaa.recipe_str,
        defaults=topoaa.params,
        mol_params=topoaa.params.pop("mol1"),
        default_params_path=None,
    )

    assert observed_inp_out == Path(protein.name).with_suffix(".inp")


def test_get_md5(topoaa, ensemble_header_w_md5, protein):
    """Test get_md5 method."""
    observed_md5_dic = topoaa.get_md5(ensemble_header_w_md5)
    expected_md5_dic = {
        1: "71098743056e0b95fbfafff690703761",
        2: "f7ab0b7c751adf44de0f25f53cfee50b",
        3: "41e028d8d28b8d97148dc5e548672142",
        4: "761cb5da81d83971c2aae2f0b857ca1e",
        5: "6c438f941cec7c6dc092c8e48e5b1c10",
    }

    assert observed_md5_dic == expected_md5_dic

    observed_md5_dic = topoaa.get_md5(protein)
    assert observed_md5_dic == {}


# The sialylation patches in carbohydrate.top all take the sialic acid as their
# "-" reference: they modify -C2 and use -C1, -O6 and -O1A, which only SIA and
# SIB have.  Their "+" reference is the acceptor, whose hydroxyl oxygen is the
# one named in the patch's `ADD BOND -C2 +O<n>` line.
SIALYL_ACCEPTOR_OXYGEN = {
    "A23": "O3",
    "A26": "O6",
    "A26S": "O6",
    "A28S": "O8",
}

OUTER_LOOP_RE = re.compile(r"^for \$id1 in id \((.*)\) loop (\w+)")
INNER_LOOP_RE = re.compile(r"^\s+for \$id2 in id \((.*)\) loop \w+")
PATCH_RE = re.compile(r'\$pres\.\$npatch="(\w+)"')


def read_bondglycans_loops():
    """Parse bondglycans.cns into one entry per outer detection loop.

    Returns a list of (loop name, $id1 selection, $id2 selection, patch names).
    """
    text = Path(RECIPE_PATH, "cns", "bondglycans.cns").read_text()
    loops = []
    for line in text.splitlines():
        outer = OUTER_LOOP_RE.match(line)
        if outer:
            loops.append((outer.group(2), outer.group(1), "", set()))
            continue
        if not loops:
            continue
        name, id1, id2, patches = loops[-1]
        inner = INNER_LOOP_RE.match(line)
        if inner:
            loops[-1] = (name, id1, inner.group(1), patches)
        patches.update(PATCH_RE.findall(line))
    return loops


def test_sialylation_loops_select_the_sialic_acid_as_first_reference():
    """Sialyl loops must put SIA/SIB in $id1, which is bound to `reference=-`.

    bondglycans.cns applies every patch with `reference=-` bound to the $id1
    residue, so a loop that selects the acceptor as $id1 hands the patch a
    residue that has none of the atoms it needs.  See issue #1711.
    """
    sialyl_loops = [
        loop
        for loop in read_bondglycans_loops()
        if loop[3] & set(SIALYL_ACCEPTOR_OXYGEN)
    ]

    assert len(sialyl_loops) == 3

    for name, id1, id2, patches in sialyl_loops:
        assert "resname SIA" in id1, name
        assert "resname SIB" in id1, name
        assert "name C2" in id1, name
        assert "resname GAL" not in id1, name
        assert "resname NGA" not in id1, name
        for patch in patches:
            assert f"name {SIALYL_ACCEPTOR_OXYGEN[patch]}" in id2, (name, patch)
