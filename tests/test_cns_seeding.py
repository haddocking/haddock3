"""CNS recipes that sample must consume the seed handed to them by Python.

``prepare_cns_input`` writes an ``evaluate ($seed=...)`` line into every ``.inp``
it builds, and modules vary it per model (``iniseed + s_ind``) so that repeated
sampling explores different trajectories. A recipe that draws random velocities
but never runs ``set seed`` silently ignores that value: CNS falls back to its
own default seed and every replica retraces the same trajectory, so
``sampling_factor > 1`` costs CPU without producing new conformations.

The check is per module rather than per file, because a module's stochastic step
and its ``set seed`` usually live in different files: ``flexref`` anneals in
``sa_ltad_*.cns`` but seeds in ``flexref.cns``, and ``rigidbody`` draws its
rotations in ``get_random_rotation.cns``. ``cgtoaa``, which inlined its MD, is
the exception rather than the rule.
"""

import re

import pytest

from haddock.modules import modules_category, modules_folder


# Velocity assignment, molecular dynamics and CNS's own random() are the
# operations that consume CNS's random number generator.
STOCHASTIC = re.compile(
    r"\bmaxwell\b|^\s*dynamics\b|\brandom\s*\(", re.IGNORECASE | re.MULTILINE
)
SET_SEED = re.compile(r"^\s*set\s+seed\b", re.IGNORECASE | re.MULTILINE)

# Every module that ships CNS recipes, taken from the registry the workflow
# engine itself validates against rather than from a filesystem glob.
CNS_MODULES = sorted(
    name
    for name, category in modules_category.items()
    if (modules_folder / category / name / "cns").is_dir()
)


def _recipes(module_name):
    """All CNS files shipped by a module, entry point and includes alike."""
    category = modules_category[module_name]
    return sorted((modules_folder / category / module_name / "cns").glob("*.cns"))


def test_cns_modules_are_discovered():
    """Guard the discovery, so the checks below cannot pass on an empty set."""
    assert {"cgtoaa", "flexref", "rigidbody", "mdref", "mdscoring"} <= set(CNS_MODULES)
    assert all(_recipes(name) for name in CNS_MODULES)


@pytest.mark.parametrize("module_name", CNS_MODULES)
def test_sampling_modules_apply_the_seed(module_name):
    """A module that draws random numbers must initialise the generator."""
    recipes = _recipes(module_name)
    sampling = [f for f in recipes if STOCHASTIC.search(f.read_text())]
    if not sampling:
        pytest.skip(f"{module_name} has no stochastic step")

    assert any(SET_SEED.search(f.read_text()) for f in recipes), (
        f"{module_name} samples in "
        f"{', '.join(f.name for f in sampling)} but no recipe runs 'set seed', "
        "so the seed prepared by Python is ignored and every replica produces "
        "the same model"
    )
