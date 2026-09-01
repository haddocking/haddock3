"""Axis 0 -- determinism of the computation itself.

Not a caching test.  It is the **precondition** for every MUST-HIT in the
matrix being meaningful, and the thing that makes conflicts interpretable.

Only artifacts that are stable bitwise, or stable after a known
normalisation, can be gated by checksum at all.  This study bounds what the
rest of the suite is permitted to assert, so it is run before Phase 2's
results are trusted -- and its findings are *reported*, not silently folded
into a pass.
"""

from __future__ import annotations

import os
from collections import defaultdict

import pytest

from cachesuite.harness import cacheable_artifacts, content_checksum, run_haddock3

pytestmark = pytest.mark.slow

#: How many times to repeat the computation.  Two is enough to find gross
#: nondeterminism; more is better and costs linearly.
REPEATS = int(os.environ.get("HADDOCK3_CACHE_AXIS0_REPEATS", "3"))


def _run(micro, tmp_path, name, **top):
    config = micro["config"].copy()
    config.top["run_dir"] = str(tmp_path / name)
    config.top.update(top)
    path = config.write(tmp_path / f"{name}.cfg")
    result = run_haddock3(path, cwd=tmp_path)
    assert result.ok, f"{name}: exit {result.returncode}\n{result.tail(30)}"
    return result


def _fingerprint(run_dir):
    return {
        artifact.relative: content_checksum(run_dir / artifact.relative)
        for artifact in cacheable_artifacts(run_dir)
    }


def test_axis0_1_repeated_runs_are_bitwise_identical(micro, tmp_path, record_property):
    """0.1 / 0.4 -- run the same jobs N times in fresh directories.

    Reported per artifact rather than as one verdict, because the useful
    answer is *which* artifacts are stable: those are the ones the rest of the
    suite may gate on.
    """
    fingerprints = [
        _fingerprint(_run(micro, tmp_path, f"determinism-{index}").run_dir)
        for index in range(REPEATS)
    ]
    reference = fingerprints[0]
    volatile = sorted(
        name
        for name in reference
        if any(other.get(name) != reference[name] for other in fingerprints[1:])
    )
    record_property("artifacts", len(reference))
    record_property("volatile", volatile)
    assert not volatile, (
        f"these artifacts are not bitwise reproducible across {REPEATS} runs "
        f"of the same configuration: {volatile}.\n"
        "This is a property of CNS and this platform, not of the cache. Until "
        "it holds, a checksum over these artifacts cannot be used to verify "
        "anything, and every MUST-HIT in the suite that touches them is "
        "unfounded."
    )


def test_axis0_2_concurrency_does_not_change_the_result(
    micro, tmp_path, record_property
):
    """0.2 -- serial versus saturated.

    If concurrency changes the bytes, then the number of cores is an
    undeclared input, and every cached result depends on the machine that
    produced it.
    """
    serial = _fingerprint(_run(micro, tmp_path, "serial", ncores=1).run_dir)
    saturated = _fingerprint(
        _run(micro, tmp_path, "saturated", ncores=max(2, os.cpu_count() or 2)).run_dir
    )
    divergent = sorted(
        name for name, value in serial.items() if saturated.get(name) != value
    )
    record_property("divergent", divergent)
    assert not divergent, (
        f"these artifacts depend on how many jobs ran at once: {divergent}. "
        "Concurrency is not in any read-set, so this is an Axis 8 defect "
        "reported by an Axis 0 experiment."
    )


def test_axis0_5_volatile_fields_are_enumerated(micro, tmp_path, record_property):
    """0.5 -- are the volatile fields exhaustively known, or only the ones
    seen so far?

    Reports *where* the differences fall when there are any, so a new volatile
    field is named rather than merely counted.
    """
    first = _run(micro, tmp_path, "fields-a").run_dir
    second = _run(micro, tmp_path, "fields-b").run_dir
    differing_lines: dict[str, list[str]] = defaultdict(list)
    for artifact in cacheable_artifacts(first):
        left = (first / artifact.relative).read_text(errors="replace").splitlines()
        right_path = second / artifact.relative
        if not right_path.exists():
            continue
        right = right_path.read_text(errors="replace").splitlines()
        for a, b in zip(left, right):
            if a != b:
                differing_lines[artifact.relative].append(f"{a!r} != {b!r}")
    record_property("differing", {k: v[:3] for k, v in differing_lines.items()})
    assert not differing_lines, (
        "unaccounted volatile content:\n"
        + "\n".join(
            f"  {name}: {examples[0]}" for name, examples in differing_lines.items()
        )
    )
