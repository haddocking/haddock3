"""What is deliberately *not* cached.

The declared boundary of what is cacheable at all -- analysis modules and
non-local execution modes are outside it.  Confirming they are **excluded**
rather than silently mishandled is a different assertion from "not tested",
and it is the one worth making: a module that is quietly served from a cache
it was never meant to use fails in exactly the way nothing else here would
catch.
"""

from __future__ import annotations

import pytest

from cachesuite.config import CACHEABLE_MODULES
from cachesuite.harness import (
    link_source,
    run_haddock3,
    step_folders,
)

pytestmark = pytest.mark.phase2


def test_analysis_modules_are_never_served_from_a_cache(
    corpus, corpus_dir, micro_config, micro, tmp_path
):
    """Analysis output must be recomputed, never linked in.

    Analysis modules are cheap, and their output describes the run they belong
    to. Serving a caprieval table from another run would attach one run's
    numbers to another run's models.
    """
    config = micro["config"].copy()
    config.top["run_dir"] = str(tmp_path / "with-analysis")
    from cachesuite.config import Step

    config.steps.append(Step("caprieval", {}))
    path = config.write(tmp_path / "with-analysis.cfg")
    source = micro["sources"]["base"]
    result = run_haddock3(
        path, cache=[source], env={"HADDOCK_CACHE_HARDLINK": "1"}, cwd=tmp_path
    )
    assert result.ok, result.tail()

    leaked = []
    for _index, module, folder in step_folders(result.run_dir):
        if module in CACHEABLE_MODULES:
            continue
        for entry in folder.rglob("*"):
            if entry.is_file() and link_source(entry, [source]) is not None:
                leaked.append(str(entry.relative_to(result.run_dir)))
    assert not leaked, (
        "these files belong to modules outside the cacheable boundary but "
        f"were served from a cache source: {leaked}"
    )


def test_non_local_execution_is_refused_rather_than_ignored(micro_config, micro, tmp_path):
    """``--cache`` outside local mode must be a defined error.

    Silently ignoring the flag would be the worst of both worlds: the user
    believes they are reusing results and pays full price, and no message says
    otherwise.
    """
    # batch mode is itself only valid with debug = true; without it the run
    # fails for that reason and the test would prove nothing about --cache.
    config = micro_config("batch-mode", mode="batch", debug=True)
    result = run_haddock3(
        config,
        cache=[micro["sources"]["base"]],
        env={"HADDOCK_CACHE_HARDLINK": "1"},
        cwd=config.parent,
    )
    assert not result.ok, (
        "--cache was accepted alongside a non-local execution mode; it must "
        "be refused with a defined error instead"
    )
    assert "cache" in result.stdout.lower(), (
        f"the error must say what was refused; output tail:\n{result.tail()}"
    )


@pytest.mark.timing
def test_cache_resolution_costs_milliseconds_per_job(corpus, record_property):
    """``t_hit`` is a requirement on the implementation, not a test margin.

    A hundred-job all-hit run must do its entire cache resolution in under a
    second, which means checksumming has to be fast or amortised.  Relaxing
    this figure 50x to accommodate slow storage does not loosen a test, it
    deletes a requirement.

    Measured as a *marginal* cost rather than a total, because a total is
    dominated by process startup and by the workflow's uncacheable remainder.
    ``tiny`` and ``tiny-wide`` are the same workflow over the same inputs and
    differ only in how many sampling jobs they schedule, so subtracting their
    all-hit reruns cancels everything except the extra jobs.  Both figures are
    measured by Phase 1 while it builds the corpus; this test does no work of
    its own beyond reading them.
    """
    from cachesuite import thresholds

    narrow, wide = corpus.require("tiny"), corpus.require("tiny-wide")
    extra_jobs = len(wide.outputs) - len(narrow.outputs)
    assert extra_jobs > 0, (
        "tiny-wide must schedule more jobs than tiny for this measurement to "
        "mean anything; the corpus is misbuilt"
    )
    marginal = (wide.allhit - narrow.allhit) / extra_jobs
    record_property("marginal_ms_per_hit", round(marginal * 1000, 2))
    record_property("extra_jobs", extra_jobs)
    budget = thresholds.T_HIT * thresholds.HIT_SCALE
    assert marginal <= budget, (
        f"resolving one cache hit costs {marginal * 1000:.1f} ms, against a "
        f"budget of {budget * 1000:.1f} ms. A hundred-job all-hit run would "
        f"spend {marginal * 100:.1f} s deciding it had nothing to do. Either "
        "the implementation is checksumming more than it needs to, or this "
        "machine is not the documented assumption (SSD, unloaded) -- see "
        "cachesuite/thresholds.py."
    )
