"""Fixtures and markers for the CNS caching contract-compliance suite.

Run order is a strict sequence, not a grouping::

    Phase 0  ->  Phase 1  ->  Phase 2

Phase 0 validates the *instrument*.  If it fails, the rest of the suite is
meaningless and must not run -- so Phase 2 aborts rather than skips.  Phase 1
is a separate script (``build_corpus.py``) whose output Phase 2 consumes.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cachesuite import damage  # noqa: E402
from cachesuite.config import Config, Step  # noqa: E402
from cachesuite.corpus import (  # noqa: E402
    PHASE0_NAME,
    CorpusNotBuilt,
    Manifest,
    assert_same_filesystem,
    corpus_root,
    freeze,
    work_root,
)
from cachesuite.generate import find_other_filesystem  # noqa: E402
from cachesuite.harness import forget_inode_index, run_haddock3  # noqa: E402
from cachesuite.systems import SYSTEMS  # noqa: E402


def pytest_configure(config):
    for marker, description in (
        ("phase0", "validates the measuring instrument; must pass before anything else"),
        ("phase2", "the suite proper; requires a corpus built by Phase 1"),
        ("timing", "asserts a Gate 2 wall-clock budget; machine-dependent"),
        ("slow", "runs CNS to completion rather than killing it at a floor"),
        ("axis", "axis(name): which taxonomy axis a case belongs to"),
    ):
        config.addinivalue_line("markers", f"{marker}: {description}")


def pytest_addoption(parser):
    group = parser.getgroup("caching suite")
    group.addoption(
        "--ignore-phase0",
        action="store_true",
        default=bool(os.environ.get("HADDOCK3_CACHE_IGNORE_PHASE0")),
        help=(
            "run Phase 2 even though Phase 0 has not passed. The instrument is "
            "then unvalidated and every Phase 2 verdict is provisional; use it "
            "only to observe how far an in-progress implementation gets."
        ),
    )
    group.addoption(
        "--no-timing",
        action="store_true",
        default=bool(os.environ.get("HADDOCK3_CACHE_NO_TIMING")),
        help=(
            "skip the Gate 2 wall-clock assertions. The correctness half of "
            "the suite -- including every MUST-MISS assertion -- is "
            "machine-independent and stays in."
        ),
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--no-timing"):
        skip = pytest.mark.skip(reason="--no-timing: Gate 2 assertions disabled")
        for item in items:
            if "timing" in item.keywords:
                item.add_marker(skip)


# --------------------------------------------------------------------------
# Phase 0 bookkeeping
# --------------------------------------------------------------------------


@pytest.fixture(scope="session")
def phase0_record() -> Path:
    """Where Phase 0's verdict is written, for Phase 2 to gate on."""
    root = corpus_root()
    root.mkdir(parents=True, exist_ok=True)
    return root / PHASE0_NAME


def pytest_sessionfinish(session, exitstatus):
    """Record whether Phase 0 passed, so a later Phase 2 session can gate."""
    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    if reporter is None:
        return
    outcomes: dict[str, str] = {}
    for outcome in ("passed", "failed", "error", "skipped"):
        for report in reporter.stats.get(outcome, []):
            nodeid = getattr(report, "nodeid", "")
            when = getattr(report, "when", "call")
            if "test_phase0" not in nodeid or when != "call":
                continue
            outcomes[nodeid] = outcome
    if not outcomes:
        return
    root = corpus_root()
    root.mkdir(parents=True, exist_ok=True)
    (root / PHASE0_NAME).write_text(
        json.dumps(
            {
                "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "outcomes": outcomes,
                "passed": all(value == "passed" for value in outcomes.values()),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )


@pytest.fixture(scope="session")
def phase0_passed(request) -> bool:
    root = corpus_root()
    path = root / PHASE0_NAME
    if not path.is_file():
        return False
    try:
        return bool(json.loads(path.read_text(encoding="utf-8")).get("passed"))
    except (ValueError, OSError):
        return False


# --------------------------------------------------------------------------
# Shared environment facts
# --------------------------------------------------------------------------


@pytest.fixture(scope="session")
def other_filesystem(tmp_path_factory) -> Path | None:
    """A writable directory on a different filesystem, or ``None``.

    Skipping the cases that need this leaves the cross-filesystem cases
    **unvalidated**, which is a coverage hole, not a neutral skip.
    """
    return find_other_filesystem(tmp_path_factory.getbasetemp())


@pytest.fixture(scope="session")
def corpus() -> Manifest:
    """The Phase 1 corpus.  Fails fast and legibly when it was never built."""
    try:
        return Manifest.read(corpus_root())
    except CorpusNotBuilt as error:
        pytest.fail(str(error), pytrace=False)


@pytest.fixture(scope="session")
def corpus_dir() -> Path:
    return corpus_root()


@pytest.fixture(scope="session")
def work_dir() -> Path:
    """A working area on the corpus's own filesystem.

    Phase 2 deliberately does not use pytest's ``tmp_path``: it lives under
    ``/tmp``, which is usually a different filesystem, and a hardlink cannot
    cross one. Every hit would arrive as a copy and be reported as a miss.
    """
    root = work_root()
    corpus = corpus_root()
    try:
        assert_same_filesystem(root, corpus)
    except RuntimeError as error:
        pytest.fail(str(error), pytrace=False)
    (root / ".gitignore").write_text("*\n", encoding="utf-8")
    # One subdirectory per pytest process, so two suites running at once
    # cannot delete each other's case directories mid-run.
    work = root / f"session-{os.getpid()}"
    work.mkdir(parents=True, exist_ok=True)
    return work


@pytest.fixture
def case_dir(work_dir, request) -> Path:
    """A private directory for one case, beside the corpus."""
    name = request.node.callspec.id if hasattr(request.node, "callspec") else request.node.name
    directory = work_dir / _safe(name)
    if directory.exists():
        shutil.rmtree(directory, onexc=_force_remove)
    directory.mkdir(parents=True)
    return directory


def _safe(name: str) -> str:
    return "".join(character if character.isalnum() or character in "-._" else "_"
                   for character in name)


def _force_remove(function, path, _excinfo):
    os.chmod(path, 0o700)
    function(path)


@pytest.fixture(scope="session", autouse=True)
def _fresh_inode_index():
    forget_inode_index()
    yield
    forget_inode_index()


# --------------------------------------------------------------------------
# The Phase 0 micro fixture -- built here, because Phase 0 runs before Phase 1
# --------------------------------------------------------------------------


@pytest.fixture(scope="session")
def micro(tmp_path_factory, other_filesystem) -> dict:
    """A minimal completed run plus the derived sources Phase 0 needs.

    Phase 0 cannot use the Phase 1 corpus: it runs first, precisely so that a
    corpus built with a broken instrument is never trusted.
    """
    root = tmp_path_factory.mktemp("micro")
    data_dir = root / "data"
    data_dir.mkdir()
    system = SYSTEMS["glycan"]
    data = {}
    for name in system.files:
        destination = data_dir / name
        shutil.copy2(system.source / name, destination)
        data[name] = destination

    base = root / "source" / "base"
    config = Config(
        top={
            "run_dir": str(base),
            "mode": "local",
            "ncores": 4,
            "clean": False,
            "postprocess": False,
            "gen_archive": False,
            "molecules": [str(data["1LMQ_r_u.pdb"]), str(data["1LMQ_l_u.pdb"])],
        },
        steps=[
            Step("topoaa"),
            Step(
                "rigidbody",
                {
                    "tolerance": 20,
                    "ambig_fname": str(data["ambig.tbl"]),
                    "sampling": 2,
                    "w_vdw": 1,
                },
            ),
        ],
    )
    config_path = config.write(root / "configs" / "micro.cfg")
    result = run_haddock3(config_path, cwd=root)
    if not result.ok:
        pytest.fail(
            "Phase 0 could not build its micro fixture; the suite cannot "
            f"measure anything.\n{result.tail(40)}",
            pytrace=False,
        )

    derived = {"base": base}

    compressed = damage.copy_run(base, root / "source" / "compressed")
    damage.compress_run(compressed)
    derived["compressed"] = compressed

    # A disjoint pair, made by deleting outputs rather than records: neither
    # source can serve the whole run, and between them they cover it.
    topo_only = damage.copy_run(base, root / "source" / "topo-only")
    damage.strip_module_artifacts(topo_only, "rigidbody")
    derived["topo-only"] = topo_only

    rigid_only = damage.copy_run(base, root / "source" / "rigid-only")
    damage.strip_module_artifacts(rigid_only, "topoaa")
    derived["rigid-only"] = rigid_only

    readonly = damage.copy_run(base, root / "source" / "readonly")
    freeze(readonly)
    derived["readonly"] = readonly

    if other_filesystem is not None:
        elsewhere = other_filesystem / "haddock3-cache-phase0"
        if elsewhere.exists():
            shutil.rmtree(elsewhere, ignore_errors=True)
        derived["crossfs"] = damage.copy_run(base, elsewhere)

    forget_inode_index()
    return {
        "root": root,
        "config": config,
        "config_path": config_path,
        "data": data,
        "sources": derived,
        "duration": result.duration,
    }


@pytest.fixture
def micro_config(micro, tmp_path):
    """A factory for B configs of the micro workflow, run out of ``tmp_path``."""

    def make(name: str = "b", **top) -> Path:
        config = micro["config"].copy()
        config.top["run_dir"] = str(tmp_path / name)
        config.top.update(top)
        return config.write(tmp_path / f"{name}.cfg")

    return make
