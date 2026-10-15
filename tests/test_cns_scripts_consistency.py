"""
Test that CNS scripts shared by several modules stay consistent.

Several modules ship their own copy of the same CNS script (e.g.
``contactairs.cns`` in flexref, emref, mdref and mdscoring). A fix applied to
one copy must be applied to all of them (see issue #1729), so every copy of a
script must be identical, unless the script is listed in
``INTENDED_VARIANTS``.

``INTENDED_VARIANTS`` records, for scripts whose copies are meant to differ,
which modules share each version. A copy that drifts from its own group, or a
variant that becomes identical to another, therefore also fails the test.
"""

from collections import defaultdict

import pytest

from haddock.modules import modules_folder


INTENDED_VARIANTS = {
    "calc_free-ene.cns": [
        {"cgtoaa", "emref", "flexref", "mdref"},
        {"emscoring", "mdscoring"},
    ],
    "flex_segment_back.cns": [
        {"cgtoaa", "emref", "flexref", "mdref"},
        {"mdscoring"},
    ],
    "flex_segment_side.cns": [
        {"flexref", "mdref"},
        {"mdscoring"},
    ],
    "generate-topology.cns": [
        {"topoaa"},
        {"topocg"},
    ],
    "print_coorheader.cns": [
        {"cgtoaa"},
        {"emref", "flexref", "mdref", "rigidbody"},
        {"emscoring", "mdscoring"},
    ],
    "read_param.cns": [
        {"cgtoaa"},
        {"emref", "flexref", "rigidbody"},
        {"emscoring"},
        {"mdref", "mdscoring"},
    ],
    "scale_inter_final.cns": [
        {"cgtoaa", "emref", "flexref", "mdref", "rigidbody"},
        {"emscoring", "mdscoring"},
    ],
    "setflags.cns": [
        {"cgtoaa", "emref", "mdref", "mdscoring"},
        {"flexref"},
    ],
}


def _shared_cns_scripts() -> dict[str, dict[str, str]]:
    """Return {script name: {module name: content}} for scripts in 2+ modules."""
    scripts: dict[str, dict[str, str]] = defaultdict(dict)
    for path in sorted(modules_folder.glob("*/*/cns/*.cns")):
        module = path.parents[1].name
        scripts[path.name][module] = path.read_text()
    return {name: copies for name, copies in scripts.items() if len(copies) > 1}


def _versions(copies: dict[str, str]) -> list[set[str]]:
    """Group the modules sharing identical content."""
    groups: dict[str, set[str]] = defaultdict(set)
    for module, content in copies.items():
        groups[content].add(module)
    return sorted(groups.values(), key=sorted)


SHARED_SCRIPTS = _shared_cns_scripts()


@pytest.mark.parametrize(
    "script",
    sorted(set(SHARED_SCRIPTS) - set(INTENDED_VARIANTS)),
)
def test_shared_cns_script_identical(script):
    """All copies of a shared CNS script are identical."""
    versions = _versions(SHARED_SCRIPTS[script])
    assert len(versions) == 1, (
        f"{script} differs between modules: {[sorted(v) for v in versions]}. "
        "Apply the change to all copies, or add the script to INTENDED_VARIANTS."
    )


@pytest.mark.parametrize("script", sorted(INTENDED_VARIANTS))
def test_cns_script_intended_variants(script):
    """Scripts with intended differences keep their expected versions."""
    expected = sorted(INTENDED_VARIANTS[script], key=sorted)
    assert _versions(SHARED_SCRIPTS.get(script, {})) == expected
