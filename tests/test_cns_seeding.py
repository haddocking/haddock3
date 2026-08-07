"""CNS recipes that sample must consume the seed handed to them by Python.

``prepare_cns_input`` writes an ``evaluate ($seed=...)`` line into every ``.inp``
it builds, and modules vary it per model (``iniseed + s_ind``) so that repeated
sampling explores different trajectories. A recipe that draws random velocities
but never runs ``set seed`` silently ignores that value: CNS falls back to its
own default seed and every replica retraces the same trajectory, so
``sampling_factor > 1`` costs CPU without producing new conformations.
"""

import re
from pathlib import Path

import pytest

from haddock import modules as modules_pkg


MODULES_PATH = Path(modules_pkg.__file__).parent

# Velocity assignment and molecular dynamics are the operations that consume
# CNS's random number generator.
STOCHASTIC = re.compile(r"\bmaxwell\b|^\s*dynamics\b", re.IGNORECASE | re.MULTILINE)
SET_SEED = re.compile(r"^\s*set\s+seed\b", re.IGNORECASE | re.MULTILINE)


def _main_recipes():
    """Yield each module's top-level CNS recipe (``<module>/cns/<module>.cns``)."""
    for cns_file in sorted(MODULES_PATH.glob("*/*/cns/*.cns")):
        module_dir = cns_file.parent.parent
        if cns_file.stem == module_dir.name:
            yield pytest.param(cns_file, id=module_dir.name)


def test_main_recipes_are_discovered():
    """Guard the glob itself, so the checks below cannot silently pass on zero."""
    found = list(_main_recipes())
    assert len(found) >= 5
    assert "cgtoaa" in {p.id for p in found}


@pytest.mark.parametrize("recipe", _main_recipes())
def test_sampling_recipes_apply_the_seed(recipe):
    """A recipe that draws random velocities must initialise the generator."""
    text = recipe.read_text()
    if not STOCHASTIC.search(text):
        pytest.skip(f"{recipe.name} has no stochastic step")

    assert SET_SEED.search(text), (
        f"{recipe.name} performs stochastic sampling but never runs "
        "'set seed', so the seed prepared by Python is ignored and every "
        "replica produces the same model"
    )
