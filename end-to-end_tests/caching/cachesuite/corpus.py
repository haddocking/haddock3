"""Locate and describe the ``OLD-RUN-DIR`` corpus.

The corpus is a **prerequisite build artifact**, in the same way a compiled
binary is: regenerable, not committed, and required before anything downstream
will work.  A checkout of the repository cannot run Phase 2.

Corpus hygiene: every fixture run directory is made read-only after it is
built, and no test ever points a run's ``run_dir`` inside one.  A test that
writes into a fixture invalidates every inode assertion that follows it.
"""

from __future__ import annotations

import json
import os
import stat
from dataclasses import asdict, dataclass, field
from pathlib import Path

SUITE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORPUS = SUITE_ROOT / "corpus"

MANIFEST_NAME = "CORPUS.json"
PHASE0_NAME = "PHASE0.json"


class CorpusNotBuilt(RuntimeError):
    """Phase 2 was reached without Phase 1 having been run."""


def corpus_root() -> Path:
    """Where the corpus lives.  ``HADDOCK3_CACHE_CORPUS`` overrides."""
    return Path(os.environ.get("HADDOCK3_CACHE_CORPUS", DEFAULT_CORPUS)).resolve()


def work_root() -> Path:
    """Where Phase 2 runs its cases.  ``HADDOCK3_CACHE_WORK`` overrides.

    **This must be on the same filesystem as the corpus.**  A hardlink cannot
    cross a filesystem boundary, so a case whose run directory lands on
    another one is served copies instead of links, Gate 1 goes blind, and
    every MUST-HIT in the case is reported as a false miss.  ``pytest``'s own
    ``tmp_path`` is under ``/tmp``, which on most machines is exactly such
    another filesystem -- which is why the suite does not use it here.
    """
    default = corpus_root().parent / "work"
    return Path(os.environ.get("HADDOCK3_CACHE_WORK", default)).resolve()


def assert_same_filesystem(work: Path, corpus: Path) -> None:
    """Fail loudly rather than reporting every hit as a miss."""
    work.mkdir(parents=True, exist_ok=True)
    if work.stat().st_dev == corpus.stat().st_dev:
        return
    raise RuntimeError(
        f"the Phase 2 working directory ({work}) and the corpus ({corpus}) "
        "are on different filesystems. Hardlinks cannot cross that boundary, "
        "so every reused result would be copied instead of linked and the "
        "whole suite would report false misses.\n"
        "Set HADDOCK3_CACHE_WORK, or HADDOCK3_CACHE_CORPUS, so that both sit "
        "on one filesystem."
    )


@dataclass
class Fixture:
    """One cache source in the corpus."""

    name: str
    kind: str  # "base" | "damaged" | "interrupted"
    path: str
    system: str = ""
    purpose: str = ""
    config: str = ""
    #: Measured wall time of the original build, seconds.
    duration: float = 0.0
    #: Measured wall time of an all-hit rerun of the same config against this
    #: fixture.  This is ``overhead(config)`` -- the only measured quantity in
    #: the timing model.  It is also the ``t_allhit`` that guards
    #: timeout-floor cases against a vacuous pass.
    allhit: float = 0.0
    outputs: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def run_dir(self, root: Path | None = None) -> Path:
        return (root or corpus_root()) / self.path


@dataclass
class Manifest:
    """Everything Phase 1 produced."""

    fixtures: dict[str, Fixture] = field(default_factory=dict)
    generated_at: str = ""
    haddock3: str = ""
    notes: list[str] = field(default_factory=list)

    def write(self, root: Path) -> None:
        payload = {
            "generated_at": self.generated_at,
            "haddock3": self.haddock3,
            "notes": self.notes,
            "fixtures": {name: asdict(f) for name, f in self.fixtures.items()},
        }
        (root / MANIFEST_NAME).write_text(
            json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
        )

    @classmethod
    def read(cls, root: Path) -> "Manifest":
        path = root / MANIFEST_NAME
        if not path.is_file():
            raise CorpusNotBuilt(
                f"corpus not built: {path} is missing.\n"
                "Phase 1 must be run before the suite proper:\n"
                "    python end-to-end_tests/caching/build_corpus.py"
            )
        payload = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            fixtures={
                name: Fixture(**data) for name, data in payload["fixtures"].items()
            },
            generated_at=payload.get("generated_at", ""),
            haddock3=payload.get("haddock3", ""),
            notes=payload.get("notes", []),
        )

    def require(self, name: str) -> Fixture:
        try:
            return self.fixtures[name]
        except KeyError:
            raise CorpusNotBuilt(
                f"corpus fixture {name!r} is missing; regenerate the corpus"
            ) from None


def freeze(path: Path) -> None:
    """Make a fixture read-only.  Directories keep ``x`` so they stay walkable."""
    for entry in sorted(path.rglob("*"), reverse=True):
        try:
            mode = entry.stat().st_mode
        except OSError:
            continue
        if entry.is_dir():
            os.chmod(entry, mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
        else:
            os.chmod(entry, mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
    os.chmod(path, path.stat().st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))


def thaw(path: Path) -> None:
    """Undo :func:`freeze` so a fixture can be replaced or deleted."""
    os.chmod(path, path.stat().st_mode | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRUSR)
    for entry in path.rglob("*"):
        try:
            os.chmod(entry, entry.stat().st_mode | stat.S_IWUSR)
        except OSError:
            continue
