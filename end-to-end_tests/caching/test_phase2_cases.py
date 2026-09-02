"""Phase 2 -- the suite proper.

Every case in ``cases/*.yaml`` is run here, through the same three gates.

Gate 1, inode identity, is primary: exact, per-job, and machine-independent.
In ``mode: complete`` the mapping is asserted **two-sided** -- a declared hit
that is missing or unlinked is a false miss (costs time), and a declared miss
that is linked is the catastrophic failure (costs correctness).

Gate 2, wall clock, is aggregate and coarse, and catches the one thing Gate 1
cannot: a "hit" that ran CNS anyway and then discarded the result.  Its
resolution degrades as the expected miss count rises, which is the second
reason Gate 1 is primary.

Gate 2 inverted, the timeout floor, is how MUST-MISS probes on expensive
modules are made affordable: run with a short ``T`` and assert the process was
still **alive** at ``T`` and had to be killed.  A self-exit before ``T`` is a
failure whatever its exit code.
"""

from __future__ import annotations

import fnmatch

import pytest

from cachesuite import thresholds
from cachesuite.cases import execute, load_cases, prepare
from cachesuite.expectations import (
    ExpectationError,
    assert_cacheable_modules_present,
    resolve,
    summarize,
)
from cachesuite.harness import (
    COMPRESSED_SUFFIXES,
    materialize_run_dir,
    cacheable_artifacts,
    content_checksum,
    exists_any,
    forget_inode_index,
    link_source,
    run_haddock3,
)

pytestmark = pytest.mark.phase2

CASES = load_cases()


def _identifier(case):
    return case.case


@pytest.fixture(scope="session", autouse=True)
def _instrument_is_validated(request, phase0_passed):
    """Phase 0 gates Phase 2.  Aborting beats skipping.

    A corpus built and measured with an unvalidated instrument produces
    verdicts that look authoritative and are not.  ``--ignore-phase0`` exists
    only so an in-progress implementation can be observed; it says so loudly.
    """
    if phase0_passed:
        return
    if request.config.getoption("--ignore-phase0"):
        rule = "=" * 72
        warning = (
            f"\n{rule}\n"
            "Phase 0 has not passed. The instrument is UNVALIDATED: "
            "HADDOCK_CACHE_HARDLINK may not be honoured, in which case Gate 1 "
            "is not a reliable hit signal and every verdict below is "
            f"provisional.\n{rule}"
        )
        request.config.pluginmanager.get_plugin("terminalreporter").write_line(
            warning, yellow=True
        )
        return
    pytest.fail(
        "Phase 0 has not passed, so the measuring instrument is unvalidated "
        "and Phase 2 results would be meaningless.\n"
        "Run it first:\n"
        "    pytest end-to-end_tests/caching/test_phase0_hardlink.py\n"
        "To observe Phase 2 anyway on an in-progress implementation, pass "
        "--ignore-phase0.",
        pytrace=False,
    )


@pytest.fixture(scope="session")
def timing_enabled(request) -> bool:
    """Whether Gate 2's wall-clock assertions apply.

    The correctness half of the suite is machine-independent -- Gate 1 is pure
    filesystem identity, valid on any hardware at any speed -- so the timing
    assertions can be dropped without weakening a single MUST-MISS.  When
    reading results from unqualified hardware: an inode failure is real
    wherever it appears, a timing failure there may be nothing.
    """
    return not request.config.getoption("--no-timing")


@pytest.mark.parametrize("case", CASES, ids=_identifier)
def test_case(
    case, corpus, corpus_dir, case_dir, other_filesystem, record_property, timing_enabled
):
    if case.skip:
        pytest.skip(f"{case.taxonomy}: {case.skip}")
    if case.needs_other_filesystem and other_filesystem is None:
        pytest.skip(
            f"COVERAGE HOLE: {case.taxonomy} needs a second filesystem and none "
            "was found; this case is unvalidated, not passing"
        )

    record_property("taxonomy", case.taxonomy)
    record_property("title", case.title)
    if case.degraded:
        record_property("degraded", case.degraded)
    if case.known_unrelated:
        record_property("known_unrelated", case.known_unrelated)

    prepared = prepare(case, corpus, corpus_dir, case_dir)
    for note in prepared.notes:
        record_property("note", note)

    if case.mode == "timeout-floor":
        _run_timeout_floor(case, prepared, record_property)
    elif case.mode == "read-set":
        _run_read_set(case, prepared, record_property, timing_enabled)
    else:
        _run_complete(case, prepared, record_property, timing_enabled)


# --------------------------------------------------------------------------
# mode: complete
# --------------------------------------------------------------------------


def _run_complete(case, prepared, record_property, timing_enabled):
    forget_inode_index()
    result = execute(prepared, timeout=_completion_timeout(case, prepared))
    record_property("duration", round(result.duration, 3))

    if result.killed:
        pytest.fail(
            f"{case.case} did not finish within "
            f"{_completion_timeout(case, prepared):.0f} s and was killed. "
            "Either far more jobs missed than the case declares, or the run "
            f"is stuck.\n{result.tail(30)}"
        )
    if case.expect_failure:
        assert not result.ok, (
            f"{case.case} was expected to fail loudly, but exited 0"
        )
        return
    assert result.ok, (
        f"{case.case}: haddock3 exited {result.returncode}\n{result.tail(40)}"
    )

    result.run_dir = materialize_run_dir(result.run_dir)
    assert_cacheable_modules_present(result.run_dir)
    try:
        expectations = resolve(result.run_dir, prepared.sources, case.expect)
    except ExpectationError as error:
        pytest.fail(f"{case.case}: {error}", pytrace=False)

    sources = list(prepared.sources.values())
    false_misses, catastrophic, wrong_entry, missing = [], [], [], []
    for expectation in expectations:
        path = result.run_dir / expectation.relative
        if not exists_any(path):
            missing.append(expectation.relative)
            continue
        actual = link_source(path, sources)
        if expectation.must_miss:
            if actual is not None:
                catastrophic.append(f"{expectation.relative} <- {actual}")
            continue
        if case.degraded:
            # Gate 1 is blind here by construction -- the result was copied
            # rather than linked, or was rewritten after the run. The
            # *catastrophic* direction above still holds, because a declared
            # miss that is nonetheless linked is unambiguous; the hit
            # direction cannot be observed at all and is left to Gate 2 and
            # Gate 3. Recorded, so the reduced power is visible in the report
            # rather than being quietly enjoyed as a pass.
            continue
        if actual is None:
            false_misses.append(
                f"{expectation.relative} (want {expectation.describe()})"
            )
        elif not _is_expected_source(actual, expectation):
            wrong_entry.append(
                f"{expectation.relative}: got {actual}, want {expectation.describe()}"
            )

    if case.degraded:
        record_property("gate1", "blind (degraded regime)")

    problems = []
    if catastrophic:
        problems.append(
            "CATASTROPHIC -- these outputs must have been recomputed but were "
            "served from the cache, so the run contains silently wrong "
            "results:\n  " + "\n  ".join(sorted(catastrophic))
        )
    if wrong_entry:
        problems.append(
            "CATASTROPHIC -- these outputs were served from the WRONG cache "
            "entry (a key collision):\n  " + "\n  ".join(sorted(wrong_entry))
        )
    if missing:
        problems.append(
            "these declared outputs are absent from the run directory:\n  "
            + "\n  ".join(sorted(missing))
        )
    if false_misses:
        problems.append(
            "FALSE MISS -- these outputs must have been reused and were "
            "recomputed instead (wasted compute, not wrong results):\n  "
            + "\n  ".join(sorted(false_misses))
        )
    if problems:
        pytest.fail(f"{case.case} ({case.taxonomy})\n\n" + "\n\n".join(problems))

    _report_renamed_hits(case, prepared, expectations, record_property)
    _gate3_content(case, prepared, result, expectations, record_property)
    _gate2_duration(case, prepared, result, expectations, record_property, timing_enabled)


def _is_expected_source(actual, expectation) -> bool:
    """Whether the file actually linked is one of the expected source entries.

    Compared modulo compression suffix.  A source stored as ``x.pdb.gz`` is
    the same cache entry as one stored as ``x.pdb`` -- the suffix is
    packaging, not identity -- and treating them as different would report
    every compressed-source case as a key collision, which is the single most
    alarming thing this suite can say.
    """

    def key(path):
        name = path.name
        for suffix in COMPRESSED_SUFFIXES:
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        return (path.parent.resolve(), name)

    return key(actual) in {key(source) for source in expectation.sources}


def _report_renamed_hits(case, prepared, expectations, record_property):
    """Count hits served from a *differently named* source file.

    This is Axis 5's evidence, and its anti-vacuity guard.  A clustering
    perturbation that happens to leave every name in place would pass the
    mapping assertion while proving nothing about rank invariance, so a case
    may require a minimum.
    """
    renamed = [
        expectation.relative
        for expectation in expectations
        if not expectation.must_miss
        and all(
            source.name != expectation.relative.rsplit("/", 1)[-1]
            for source in expectation.sources
        )
    ]
    record_property("renamed_hits", len(renamed))
    if case.require_renamed_hits:
        assert len(renamed) >= case.require_renamed_hits, (
            f"{case.case}: only {len(renamed)} of the hits came from a "
            f"differently named source file, but the case needs at least "
            f"{case.require_renamed_hits}. The perturbation did not actually "
            "reorder or rename anything, so this case proves nothing -- fix "
            "the case, do not relax the requirement."
        )


def _run_read_set(case, prepared, record_property, timing_enabled):
    """Axis 8 and Axis 10 -- is a hit here *sound*?

    A pairwise MUST-HIT/MUST-MISS assertion cannot see an undeclared
    dependency: if the perturbation influences the result but appears in no
    read-set, the key does not change and the job hits, which looks correct.
    So this mode runs B twice -- once against the cache, once with no cache at
    all -- and asserts that what the cache served is byte-identical to what a
    fresh computation under the same perturbation produces.

    If they differ, the perturbation reached the result without reaching the
    key.  That is the catastrophic failure, and it is invisible to every other
    mode in this suite.
    """
    forget_inode_index()
    cached = execute(prepared, timeout=_completion_timeout(case, prepared))
    assert cached.ok, f"{case.case} (cached): exit {cached.returncode}\n{cached.tail(30)}"
    record_property("duration", round(cached.duration, 3))

    uncached_dir = prepared.run_dir.parent / f"{prepared.run_dir.name}-uncached"
    text = prepared.config_path.read_text(encoding="utf-8").replace(
        f'run_dir = "{prepared.run_dir}"', f'run_dir = "{uncached_dir}"'
    )
    uncached_config = prepared.config_path.with_name("uncached.cfg")
    uncached_config.write_text(text, encoding="utf-8")
    reference = run_haddock3(
        uncached_config,
        cache=[],
        env=prepared.env,
        cwd=prepared.cwd,
        install=prepared.install,
        timeout=_completion_timeout(case, prepared),
    )
    assert reference.ok, (
        f"{case.case} (uncached reference): exit {reference.returncode}\n"
        f"{reference.tail(30)}"
    )

    divergent = []
    for artifact in cacheable_artifacts(cached.run_dir):
        fresh = uncached_dir / artifact.relative
        if not exists_any(fresh):
            continue
        if content_checksum(cached.run_dir / artifact.relative) != content_checksum(fresh):
            divergent.append(artifact.relative)
    assert not divergent, (
        f"{case.case} ({case.taxonomy}): CATASTROPHIC -- the cache served "
        f"results that differ from what this perturbation actually computes: "
        f"{divergent}. The perturbation reached the result without reaching "
        "the key, so the read-set is incomplete (Axis 8) and every hit this "
        "key has ever served was already suspect."
    )


def _completion_timeout(case, prepared) -> float:
    """A ceiling that stops a stuck run.  Gate 2 is the real bound.

    Anchored on how long the base run took to build *from scratch*, because a
    case is a variation on that workflow: even if every job misses, it is
    doing at most what the base run already did.  A flat half-hour ceiling
    would be correct but would let one hung case cost more wall time than the
    rest of the suite put together.
    """
    misses = case.misses or {}
    declared = sum(
        thresholds.miss_budget(module) * count for module, count in misses.items()
    )
    scale = thresholds.MISS_SCALE
    return max((prepared.base_duration * 3 + declared) * scale, 300.0 * scale)


def _gate2_duration(case, prepared, result, expectations, record_property, enabled):
    """Gate 2 -- the declared-budget upper bound.

    Recorded on pass as well as fail, so drift is visible before it becomes a
    flake.
    """
    misses, hits = summarize(expectations)
    if case.misses is not None:
        misses = dict(case.misses)
    bound = thresholds.gate2_bound(prepared.overhead, misses, hits)
    record_property("gate2_bound", round(bound, 2))
    if not enabled:
        return
    assert result.duration <= bound, (
        f"{case.case}: took {result.duration:.1f} s, budget {bound:.1f} s "
        f"(overhead {prepared.overhead:.1f} + misses {misses} + {hits} hits). "
        "Either a job that should have hit ran CNS, or this machine is slower "
        "than the documented assumption (SSD, unloaded); see thresholds.py."
    )


def _gate3_content(case, prepared, result, expectations, record_property):
    """Gate 3 -- content equality.

    Byte-identity is **not** evidence of a hit: a deterministic recomputation
    produces the same bytes.  It is used only to confirm that a *copy* rather
    than a link delivered the right bytes, which is all that is left in the
    degraded copy-mode regime.
    """
    if not case.degraded:
        return
    # In the degraded regime this is the only evidence that a hit delivered
    # the right result at all, so it is asserted here and nowhere else.
    mismatched = []
    for expectation in expectations:
        if expectation.must_miss:
            continue
        got = content_checksum(result.run_dir / expectation.relative)
        wanted = {
            content_checksum(source)
            for source in expectation.sources
            if exists_any(source)
        }
        if wanted and got not in wanted:
            mismatched.append(expectation.relative)
    assert not mismatched, (
        f"{case.case}: copy-mode delivered the wrong bytes for {mismatched}"
    )


# --------------------------------------------------------------------------
# mode: timeout-floor
# --------------------------------------------------------------------------


def _run_timeout_floor(case, prepared, record_property):
    """Gate 2 inverted.

    Run a MUST-MISS case with a deliberately short ``T`` and assert the run had
    to be killed:

    * the job genuinely missed -> CNS starts -> still running at ``T`` -> pass
    * the job wrongly hit -> the run completes in milliseconds -> fail

    Choosing ``T`` is easy rather than delicate: hits cost ~10 ms and misses
    cost seconds to minutes, so the window is three to four orders of
    magnitude wide.
    """
    floor = case.floor or thresholds.DEFAULT_TIMEOUT_FLOOR
    allhit = prepared.overhead
    record_property("floor", floor)
    record_property("t_allhit", round(allhit, 2))

    # Anti-vacuity guard, primary: if T is shorter than the time to *reach*
    # the perturbed job, the run times out before caching was ever consulted
    # and the test passes having proved nothing. The all-hit run of the
    # unperturbed config reaches every job; if it finishes well inside T, then
    # the perturbed run still being alive at T means it did something the
    # unperturbed run did not, which is CNS execution. This makes no
    # assumption about ordering or resolution strategy -- resolution may be
    # interleaved with execution.
    if allhit > floor * thresholds.ALLHIT_FLOOR_RATIO:
        pytest.fail(
            f"{case.case}: vacuity guard. The unperturbed all-hit run of this "
            f"config takes {allhit:.1f} s, which is not far enough below the "
            f"{floor:.1f} s floor for 'still alive at T' to mean 'ran CNS'. "
            "Raise floor for this case, or move it to mode: complete.",
            pytrace=False,
        )

    forget_inode_index()
    result = execute(prepared, timeout=floor)
    record_property("duration", round(result.duration, 3))

    assert result.killed, (
        f"{case.case} ({case.taxonomy}): the run exited on its own after "
        f"{result.duration:.2f} s, before the {floor:.1f} s floor. A "
        "MUST-MISS job would still have been running CNS. Either it wrongly "
        "hit -- the catastrophic failure -- or it crashed; the two look alike "
        f"from outside and need investigation to separate.\n{result.tail(30)}"
    )

    # Partial-directory inspection. In this mode, and only in this mode, Gate 1
    # is one-sided: a declared-miss output that is present and linked is still
    # a failure, but the *absence* of a declared hit is indistinguishable from
    # not-yet-reached under interleaved resolution.
    sources = list(prepared.sources.values())
    if result.run_dir.is_dir():
        linked = [
            artifact.relative
            for artifact in cacheable_artifacts(result.run_dir)
            if _declared_miss(case, artifact.relative)
            and link_source(result.run_dir / artifact.relative, sources) is not None
        ]
        assert not linked, (
            f"{case.case}: these outputs must have been recomputed but are "
            f"already hardlinked into a cache source: {linked}"
        )

    # Earlier-step witness, reinforcement. HADDOCK3 runs steps sequentially,
    # so a MUST-HIT job in a strictly earlier step must be present and linked
    # at kill time. Within a step, job order is not guaranteed under the
    # parallel scheduler, so a witness is never a sibling.
    if case.witness and result.run_dir.is_dir():
        witnesses = [
            artifact
            for artifact in cacheable_artifacts(result.run_dir)
            if fnmatch.fnmatch(artifact.relative, case.witness)
        ]
        assert witnesses, (
            f"{case.case}: the witness {case.witness!r} matched nothing in the "
            "partial run directory, so reaching the perturbed job was never "
            "established"
        )
        unlinked = [
            artifact.relative
            for artifact in witnesses
            if link_source(result.run_dir / artifact.relative, sources) is None
        ]
        assert not unlinked, (
            f"{case.case}: the earlier-step witness was not served from the "
            f"cache: {unlinked}. Either the run never reached the perturbed "
            "job, or an unrelated MUST-HIT is failing."
        )


def _declared_miss(case, relative: str) -> bool:
    paths = case.expect.get("paths") or {}
    if relative in paths:
        return paths[relative] is None
    for rule in case.expect.get("patterns") or []:
        if fnmatch.fnmatch(relative, rule["match"]):
            return rule["expect"] == "miss"
    return case.expect.get("default", "hit") == "miss"
