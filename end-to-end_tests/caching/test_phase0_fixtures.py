"""Phase 0 -- does each fixture produce the situation its case describes?

Part of the instrument, not of the feature: nothing here runs ``haddock3``.

A case is only as good as its perturbation. A fixture that leaves the file it
claims to change byte-identical, or that changes it in a way the case does not
describe, produces a verdict that looks authoritative and is not -- and it
fails, or passes, for a reason the case never names. Those are not
hypothetical: several of them were shipped in this suite and were invisible
until an implementation disagreed with them.

So the perturbations whose *exact shape* a case's verdict depends on are
checked here, mechanically, before any of them is used to judge anything.
"""

from pathlib import Path

import pytest

from cachesuite.cases import load_cases
from cachesuite.perturb import (
    Context,
    PerturbationError,
    _edit_models,
    _install_replace,
)
from cachesuite.systems import REPO_ROOT

pytestmark = pytest.mark.phase0


def _ensemble(tmp_path: Path, members: int = 4) -> Path:
    path = tmp_path / "ensemble.pdb"
    blocks = []
    for index in range(1, members + 1):
        blocks.append(
            f"MODEL     {index}\n"
            f"ATOM      1  CA  ALA A   1    {index:8.3f}   0.000   0.000\n"
            "ENDMDL\n"
        )
    path.write_text("REMARK  an ensemble\n" + "".join(blocks), encoding="utf-8")
    return path


def _members(path: Path) -> list[str]:
    import re

    return re.findall(r"MODEL.*?ENDMDL\s*\n", path.read_text(encoding="utf-8"), re.S)


def _context(tmp_path: Path) -> Context:
    return Context(config=None, data={}, data_dir=tmp_path, tmp=tmp_path)


def test_removing_a_member_takes_one_from_the_middle(tmp_path):
    """Otherwise every remaining name still matches its old content.

    Dropping the last member is not the perturbation axis6.3 describes: it
    leaves the surviving members exactly where they were, so the case would
    pass without exercising the renumbering it is named after.
    """
    path = _ensemble(tmp_path)
    before = _members(path)

    _edit_models(path, {"models": "remove"}, _context(tmp_path))

    after = _members(path)
    assert len(after) == len(before) - 1
    assert after[-1] == before[-1], "the last member must survive"
    assert after != before[:-1], "a middle member must be the one removed"


def test_adding_a_member_appends_a_byte_identical_copy(tmp_path):
    """axis6.2b's premise: the added member is not new work at all."""
    path = _ensemble(tmp_path)
    before = _members(path)

    _edit_models(path, {"models": "add"}, _context(tmp_path))

    after = _members(path)
    assert after[:-1] == before
    assert after[-1] == before[0]


def test_adding_a_distinct_member_appends_a_new_conformer(tmp_path):
    """axis6.2a's premise: the added member is a computation nobody has done."""
    path = _ensemble(tmp_path)
    before = _members(path)

    _edit_models(path, {"models": "add_distinct"}, _context(tmp_path))

    after = _members(path)
    assert after[:-1] == before
    assert after[-1] not in before


def test_reordering_moves_members_without_changing_any(tmp_path):
    """axis6.4's premise: same bytes, different index."""
    path = _ensemble(tmp_path)
    before = _members(path)

    _edit_models(path, {"models": "reorder"}, _context(tmp_path))

    after = _members(path)
    assert sorted(after) == sorted(before)
    assert after[0] == before[1] and after[1] == before[0]


def test_a_recipe_substitution_must_match_exactly_once(tmp_path):
    """A silent no-op would invert a case's meaning, so it is an error."""
    path = tmp_path / "recipe.cns"

    path.write_text("nothing to see\n", encoding="utf-8")
    with pytest.raises(PerturbationError):
        _install_replace(path, "harm = 1", "harm = 5")

    path.write_text("harm = 1\nharm = 1\n", encoding="utf-8")
    with pytest.raises(PerturbationError):
        _install_replace(path, "harm = 1", "harm = 5")

    path.write_text("harm = 1\n", encoding="utf-8")
    _install_replace(path, "harm = 1", "harm = 5")
    assert path.read_text(encoding="utf-8") == "harm = 5\n"


def test_every_declared_recipe_substitution_still_applies():
    """A case that substitutes a line must find that line exactly once.

    Read from the cases themselves rather than restated here, so the check
    cannot drift away from what is actually asserted. If a recipe is reworded,
    this fails loudly in the instrument rather than reducing a case that
    claims a substantive change to a second copy of the case that claims an
    inert one -- asserting the opposite verdict, and passing for the wrong
    reason.
    """
    checked = 0
    for case in load_cases():
        for operation in case.perturb:
            if operation.get("op") != "install":
                continue
            for relative, how in (operation.get("edit") or {}).items():
                if not isinstance(how, dict):
                    continue
                wanted, _replacement = how["replace"]
                recipe = (REPO_ROOT / "src/haddock" / relative).read_text("utf-8")
                assert recipe.count(wanted) == 1, (
                    f"{case.case} substitutes a line of {relative} that occurs "
                    f"{recipe.count(wanted)} times; the fixture no longer makes "
                    "the change the case describes"
                )
                checked += 1
    assert checked, "no case declares a recipe substitution any more"
