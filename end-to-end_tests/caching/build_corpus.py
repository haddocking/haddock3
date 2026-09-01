#!/usr/bin/env python3
"""Phase 1 entry point: build the ``OLD-RUN-DIR`` corpus.

The corpus is a prerequisite build artifact.  Phase 2 cannot run without it,
and this script is the only thing that produces it::

    python end-to-end_tests/caching/build_corpus.py
    python end-to-end_tests/caching/build_corpus.py --only tiny --only pp
    HADDOCK3_CACHE_CORPUS=/scratch/corpus python .../build_corpus.py

Expensive by design: it runs real CNS.  Everything downstream is cheap because
this is paid once.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cachesuite import generate  # noqa: E402
from cachesuite.corpus import corpus_root  # noqa: E402
from cachesuite.systems import BASE_RUNS  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="corpus root (default: $HADDOCK3_CACHE_CORPUS or ./corpus)",
    )
    parser.add_argument(
        "--only",
        action="append",
        default=[],
        choices=[spec.name for spec in BASE_RUNS],
        help="build only these base runs, and no derived fixtures; repeatable",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "skip fixtures already recorded in the manifest and still present "
            "on disk; the build writes the manifest after every fixture, so "
            "an interrupted run can be continued with this flag"
        ),
    )
    parser.add_argument(
        "--skip-interrupted",
        action="store_true",
        help="skip the Axis 9 fixtures, which are the least deterministic",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="describe the corpus without building it",
    )
    arguments = parser.parse_args(argv)

    if arguments.list:
        for spec in BASE_RUNS:
            print(f"{spec.name:<14} [{spec.cost}] {spec.system}")
            print(f"{'':<14} {spec.purpose}")
        return 0

    root = (arguments.root or corpus_root()).resolve()
    print(f"[corpus] root: {root}")
    manifest = generate.build(
        root,
        only=tuple(arguments.only) or None,
        skip_interrupted=arguments.skip_interrupted,
        resume=arguments.resume,
    )
    print(f"[corpus] {len(manifest.fixtures)} fixtures")
    for note in manifest.notes:
        print(f"[corpus] {note}")
    unusable = [f.name for f in manifest.fixtures.values() if not f.path]
    if unusable:
        print(f"[corpus] unusable fixtures: {', '.join(unusable)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
