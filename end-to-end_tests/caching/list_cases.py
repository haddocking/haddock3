#!/usr/bin/env python3
"""Print the case matrix in readable form.

The YAML in ``cases/`` is the specification.  This just renders it, so that
reviewing what is claimed does not require reading YAML::

    python end-to-end_tests/caching/list_cases.py
    python end-to-end_tests/caching/list_cases.py --markdown > CASES.md
    python end-to-end_tests/caching/list_cases.py --axis 5
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cachesuite.cases import load_cases  # noqa: E402


def _wrap(text: str, indent: str, width: int = 78) -> str:
    return textwrap.fill(
        " ".join(text.split()),
        width=width,
        initial_indent=indent,
        subsequent_indent=indent,
    )


def plain(cases) -> None:
    current = None
    for case in cases:
        if case.source_file != current:
            current = case.source_file
            print(f"\n{'=' * 78}\n{current}\n{'=' * 78}")
        marker = "SKIPPED" if case.skip else case.mode
        print(f"\n{case.taxonomy or '-':<12} {case.case}   [{marker}]")
        if case.title:
            print(_wrap(case.title, "    "))
        if case.why:
            print(_wrap(case.why, "      | "))
        if case.skip:
            print(_wrap(f"NOT RUN: {case.skip}", "      ! "))
        if case.degraded:
            print(_wrap(f"REDUCED POWER: {case.degraded}", "      ~ "))


def markdown(cases) -> None:
    print("# CNS caching: case matrix\n")
    print(f"{len(cases)} cases, "
          f"{sum(1 for case in cases if not case.skip)} of them executable.\n")
    current = None
    for case in cases:
        if case.source_file != current:
            current = case.source_file
            print(f"\n## `{current}`\n")
            print("| Taxonomy | Case | Mode | What it claims |")
            print("|---|---|---|---|")
        why = " ".join((case.title + ". " + case.why).split())
        if case.skip:
            why = "**Not run.** " + " ".join(case.skip.split())
        print(
            f"| {case.taxonomy or ''} | `{case.case}` | "
            f"{'skipped' if case.skip else case.mode} | {why} |"
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markdown", action="store_true", help="emit a table")
    parser.add_argument("--axis", default=None, help="only this axis")
    parser.add_argument(
        "--runnable", action="store_true", help="hide cases marked out of scope"
    )
    arguments = parser.parse_args(argv)

    cases = load_cases()
    if arguments.axis:
        cases = [
            case
            for case in cases
            if case.taxonomy.startswith(arguments.axis)
            or case.source_file.startswith(f"axis{arguments.axis}")
        ]
    if arguments.runnable:
        cases = [case for case in cases if not case.skip]

    (markdown if arguments.markdown else plain)(cases)
    if not arguments.markdown:
        executable = sum(1 for case in cases if not case.skip)
        print(
            f"\n{len(cases)} cases: {executable} executable, "
            f"{len(cases) - executable} recorded as out of scope."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
