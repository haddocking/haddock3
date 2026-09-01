"""Phase 0 -- validate the instrument, before anything is measured with it.

These tests do not test the cache.  They test ``HADDOCK_CACHE_HARDLINK``, the
switch that makes Gate 1 -- inode identity -- a *total* observable.  If they
fail, every verdict in the rest of the suite is unreliable, so Phase 2 aborts
rather than skips.

Why an environment variable and not a config parameter: a config parameter
would be inside the config file and would therefore have to be *tested* for
key-invariance (Axis 3).  An environment variable makes that
untestable-because-impossible -- discharge by construction, applied to the test
instrument itself.

===========  ==========================================================
Value        Meaning
===========  ==========================================================
undefined    best-effort: hardlink where possible, copy otherwise
``0``        force copy, always
``1``        force hardlink; **any** hardlinking failure is an error
===========  ==========================================================

The ``1`` case matters more than it looks.  A silent fallback to copy under
``=1`` would report every MUST-HIT in the entire suite as a false MISS.
Failing loudly is what keeps the instrument honest.
"""

from __future__ import annotations

import pytest

from cachesuite.harness import (
    cacheable_artifacts,
    content_checksum,
    exists_any,
    forget_inode_index,
    link_source,
    run_haddock3,
)

pytestmark = pytest.mark.phase0

HARDLINK = "HADDOCK_CACHE_HARDLINK"


def _run(micro_config, micro, sources, value, name="b", **top):
    env = {} if value is None else {HARDLINK: value}
    config = micro_config(name, **top)
    forget_inode_index()
    return run_haddock3(
        config,
        cache=[micro["sources"][key] for key in sources],
        env=env,
        cwd=config.parent,
    )


def _source_paths(micro, sources):
    return [micro["sources"][key] for key in sources]


def _linked(result, micro, sources):
    """``{relative: source path or None}`` for every cacheable output."""
    roots = _source_paths(micro, sources)
    return {
        artifact.relative: link_source(result.run_dir / artifact.relative, roots)
        for artifact in cacheable_artifacts(result.run_dir)
    }


def test_p0_1_same_filesystem_forced_hardlink_shares_inodes(micro, micro_config):
    """P0.1 -- same fs, ``=1``, full cache: every reused file shares an inode."""
    result = _run(micro_config, micro, ["base"], "1")
    assert result.ok, result.tail()
    links = _linked(result, micro, ["base"])
    assert links, "the run produced no cacheable outputs at all"
    unlinked = sorted(name for name, source in links.items() if source is None)
    assert not unlinked, (
        "under HADDOCK_CACHE_HARDLINK=1 with a complete cache, every output "
        f"must be a hardlink into the source; these were not: {unlinked}"
    )


def test_p0_2_forced_copy_does_not_share_inodes(micro, micro_config):
    """P0.2 -- same fs, ``=0``: inodes differ, content identical, run succeeds."""
    result = _run(micro_config, micro, ["base"], "0")
    assert result.ok, result.tail()
    links = _linked(result, micro, ["base"])
    linked = sorted(name for name, source in links.items() if source is not None)
    assert not linked, (
        "HADDOCK_CACHE_HARDLINK=0 must force a copy, but these outputs share "
        f"an inode with the source: {linked}"
    )
    source = micro["sources"]["base"]
    for artifact in cacheable_artifacts(result.run_dir):
        expected = source / artifact.relative
        if exists_any(expected):
            assert content_checksum(
                result.run_dir / artifact.relative
            ) == content_checksum(expected), (
                f"{artifact.relative} was copied with the wrong bytes"
            )


def test_p0_3_cross_filesystem_forced_hardlink_fails_loudly(
    micro, micro_config, other_filesystem
):
    """P0.3 -- cross fs, ``=1``: a defined error, never a silent copy.

    This is the load-bearing one.  A silent fallback to copy here is exactly
    what would turn every MUST-HIT in the suite into a false MISS while
    reporting success.
    """
    if other_filesystem is None:
        pytest.skip(
            "COVERAGE HOLE: no second filesystem available, so the suite's "
            "cross-filesystem cases are unvalidated"
        )
    result = _run(micro_config, micro, ["crossfs"], "1")
    assert not result.ok, (
        "HADDOCK_CACHE_HARDLINK=1 with a cache source on another filesystem "
        "must fail loudly; the run exited 0 instead, which means it fell back "
        "to copying silently"
    )


def test_p0_4_cross_filesystem_forced_copy_succeeds(
    micro, micro_config, other_filesystem
):
    """P0.4 -- cross fs, ``=0``: copies, run succeeds."""
    if other_filesystem is None:
        pytest.skip("COVERAGE HOLE: no second filesystem available")
    result = _run(micro_config, micro, ["crossfs"], "0")
    assert result.ok, result.tail()
    links = _linked(result, micro, ["crossfs"])
    assert not [name for name, source in links.items() if source is not None]


def test_p0_5_cross_filesystem_best_effort_copies(
    micro, micro_config, other_filesystem
):
    """P0.5 -- cross fs, undefined: best-effort degrades to copying."""
    if other_filesystem is None:
        pytest.skip("COVERAGE HOLE: no second filesystem available")
    result = _run(micro_config, micro, ["crossfs"], None)
    assert result.ok, result.tail()
    links = _linked(result, micro, ["crossfs"])
    assert not [name for name, source in links.items() if source is not None]


def test_p0_6_same_filesystem_best_effort_hardlinks(micro, micro_config):
    """P0.6 -- same fs, undefined: best-effort takes the hardlink."""
    result = _run(micro_config, micro, ["base"], None)
    assert result.ok, result.tail()
    links = _linked(result, micro, ["base"])
    assert any(source is not None for source in links.values()), (
        "with no HADDOCK_CACHE_HARDLINK set and a same-filesystem source, "
        "best-effort must produce at least some hardlinks"
    )


@pytest.mark.parametrize("value", ["2", "yes", "", "true", "-1"])
def test_p0_7_invalid_value_is_a_defined_error(micro, micro_config, value):
    """P0.7 -- an invalid value is rejected, not silently defaulted.

    A silent fall-through to best-effort would mean a typo in the suite's own
    environment quietly disables Gate 1 for every case that follows.
    """
    result = _run(micro_config, micro, ["base"], value, name=f"b{abs(hash(value))}")
    assert not result.ok, (
        f"HADDOCK_CACHE_HARDLINK={value!r} is not one of 0/1/unset and must be "
        "rejected with a defined error; the run exited 0 instead"
    )
    assert HARDLINK in result.stdout, (
        "the error must name HADDOCK_CACHE_HARDLINK, so the operator can see "
        f"what was rejected; output tail:\n{result.tail()}"
    )


def test_p0_8_compressed_source_either_links_the_archive_or_errors(
    micro, micro_config, record_property
):
    """P0.8 -- pin down which regime the compressed-source cases run in.

    Either a ``.pdb.gz`` source is hardlinked *as such* (Gate 1 stays live, and
    Axes 1.7 / 12.7 stay in the strong regime), or the implementation must
    decompress and ``=1`` errors (those axes move to copy mode, where Gate 1 is
    blind).  Both are defensible; a silent copy under ``=1`` is not.
    """
    result = _run(micro_config, micro, ["compressed"], "1")
    if not result.ok:
        record_property("p0_8_regime", "errors")
        pytest.skip(
            "P0.8: this implementation cannot hardlink a compressed source "
            "under =1. Axes 1.7 and 12.7 therefore run in the degraded "
            "copy-mode regime, where Gate 1 is blind."
        )
    links = _linked(result, micro, ["compressed"])
    linked = [name for name, source in links.items() if source is not None]
    record_property("p0_8_regime", "links-the-archive" if linked else "silent-copy")
    assert linked, (
        "the run succeeded under HADDOCK_CACHE_HARDLINK=1 with a compressed "
        "source but produced no hardlink into it, which means it copied "
        "silently -- the one outcome =1 forbids"
    )


def test_p0_9_read_only_source_can_still_be_linked(micro, micro_config):
    """P0.9 -- a read-only source links fine: linking needs target-dir write."""
    result = _run(micro_config, micro, ["readonly"], "1")
    assert result.ok, result.tail()
    links = _linked(result, micro, ["readonly"])
    unlinked = sorted(name for name, source in links.items() if source is None)
    assert not unlinked, (
        "a hardlink needs write permission on the *target* directory, not on "
        f"the source; these outputs were not linked: {unlinked}"
    )


def test_p0_10_multiple_sources_resolve_from_whichever_holds_the_entry(
    micro, micro_config
):
    """P0.10 -- ``--cache`` is repeatable and coverage is the union."""
    sources = ["topo-only", "rigid-only"]
    result = _run(micro_config, micro, sources, "1")
    assert result.ok, result.tail()
    links = _linked(result, micro, sources)
    unlinked = sorted(name for name, source in links.items() if source is None)
    assert not unlinked, (
        "the union of two disjoint sources covers every job; these were not "
        f"resolved from either: {unlinked}"
    )
    by_source = {
        name: source.parent.parent.name for name, source in links.items() if source
    }
    assert set(by_source.values()) == set(sources), (
        "each source must actually be used for the half it holds; resolution "
        f"came from {sorted(set(by_source.values()))}"
    )
