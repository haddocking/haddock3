"""Declared timing budgets for Gate 2.

These are **documented assumptions, not measurements**: defaults for a
reasonable machine with an SSD that is not heavily loaded.  There is no
autocalibration and no machine profile.  Anyone running elsewhere adjusts the
constants here (or the environment overrides) and knows they have done so.

``T_HIT`` is not an arbitrary margin.  It is a genuine performance requirement
on the implementation: a 100-job all-hit run must resolve in under a second.
Relaxing it 50x to accommodate slow storage does not loosen a test, it deletes
a requirement.
"""

from __future__ import annotations

import os

#: Seconds per cache miss, by module.  Generous declared ceilings, never
#: estimates, so they stay reproducible across machines without margin tuning.
T_MISS: dict[str, float] = {
    "topoaa": 10.0,
    "topocg": 10.0,
    "rigidbody": 10.0,
    "cgtoaa": 10.0,
    "flexref": 180.0,
    "emref": 180.0,
    "mdref": 180.0,
    "emscoring": 180.0,
    "mdscoring": 180.0,
}

#: Seconds per cache hit, **averaged** over the run.  Per-job jitter is fine.
T_HIT = 0.010

#: Fallback when ``overhead`` was never measured for a config.
DEFAULT_OVERHEAD = 60.0

#: ``overhead`` is the one *measured* quantity in the timing model, and it is
#: measured once.  Python startup, module import and filesystem cache state
#: vary by hundreds of milliseconds between otherwise identical runs, so a
#: bound built on a single measurement with no slack fails on jitter alone --
#: which is not a caching defect and must not be reported as one.
#:
#: These two absorb that.  They are deliberately modest: with the tiny base
#: run (overhead ~4 s, eight jobs) the bound lands near 7 s, while recomputing
#: all eight jobs would cost ~12 s.  Gate 2 therefore still catches a run that
#: executed CNS when it should have reused, which is the one thing it exists
#: for.  Widening them further would make it catch nothing.
OVERHEAD_TOLERANCE = 1.25
OVERHEAD_SLACK = 2.0

#: Default ``T`` for Gate 2 inverted.  The window between ``t_allhit`` and
#: ``t_miss`` is three to four orders of magnitude wide, so this is a choice
#: rather than a calibration.
DEFAULT_TIMEOUT_FLOOR = 8.0

#: ``t_allhit`` must be at most this fraction of the floor for a timeout-floor
#: case to be non-vacuous.
ALLHIT_FLOOR_RATIO = 0.5

#: Multiplies every ``T_MISS`` and ``overhead``.  The one knob for slow CPUs.
MISS_SCALE = float(os.environ.get("HADDOCK3_CACHE_MISS_SCALE", "1.0"))

#: Multiplies ``T_HIT``.  Separated from ``MISS_SCALE`` on purpose: raising
#: this is deleting a requirement, and should be a deliberate, visible act.
HIT_SCALE = float(os.environ.get("HADDOCK3_CACHE_HIT_SCALE", "1.0"))


def miss_budget(module: str) -> float:
    """Declared ceiling for one cache miss in ``module``."""
    return T_MISS.get(module, 180.0) * MISS_SCALE


def hit_budget(n_hits: int) -> float:
    """Declared ceiling for ``n_hits`` cache hits, in aggregate."""
    return n_hits * T_HIT * HIT_SCALE


def gate2_bound(overhead: float, misses: dict[str, int], n_hits: int) -> float:
    """``overhead + sum(t_miss) + n_hits * t_hit`` -- the Gate 2 upper bound.

    The measured overhead carries a tolerance; the declared budgets do not,
    because they are already generous ceilings rather than estimates.
    """
    total = (overhead * OVERHEAD_TOLERANCE + OVERHEAD_SLACK) * MISS_SCALE
    total += hit_budget(n_hits)
    for module, count in misses.items():
        total += miss_budget(module) * count
    return total
