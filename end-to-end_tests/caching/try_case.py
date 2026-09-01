#!/usr/bin/env python3
"""Run, read, or export any single caching test as ordinary HADDOCK3.

The test suite runs under pytest, which is convenient for running 140 cases at
once and inconvenient for understanding any one of them.  This script is the
other door.  It takes one case and gives you:

* ``--explain``  what the case does and what should happen, in plain words
* ``--run``      the case, run once, with a table of what was reused
* ``--export``   a self-contained folder with two .cfg files and a run.sh,
                 which you can run, read and edit without this suite at all

The exported folder is the point.  It contains nothing but HADDOCK3
configuration files and a shell script -- the same things you would write by
hand -- so you can change a parameter, rerun, and see what happens to the
reuse.  Nothing in it imports the test harness.

Examples::

    python try_case.py --list
    python try_case.py --list --common
    python try_case.py axis4.1-sampling-increased
    python try_case.py axis4.1-sampling-increased --run
    python try_case.py axis4.1-sampling-increased --export ~/cache-demo
"""

from __future__ import annotations

import argparse
import shutil
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cachesuite.cases import Case, load_cases, prepare  # noqa: E402
from cachesuite.corpus import CorpusNotBuilt, Manifest, corpus_root  # noqa: E402
from cachesuite.harness import cacheable_artifacts, link_source  # noqa: E402

#: The standalone checker in ``by-hand/``.  It imports nothing from this suite
#: on purpose, which is what makes an exported case runnable without it.
CHECK_REUSE = Path(__file__).resolve().parent / "by-hand" / "check_reuse.py"


def _wrap(text: str, indent: str = "  ") -> str:
    return "\n".join(
        textwrap.fill(
            paragraph.strip(),
            width=78,
            initial_indent=indent,
            subsequent_indent=indent,
        )
        for paragraph in text.strip().split("\n\n")
        if paragraph.strip()
    )


def _plain_expectation(case: Case) -> str:
    """State the case's verdict in words rather than in schema."""
    expect = case.expect or {}
    default = expect.get("default", "hit")
    lines = []
    if case.skip:
        return "This case is deliberately not run. " + " ".join(case.skip.split())
    if case.expect_failure:
        lines.append("HADDOCK3 should refuse to run and say why.")
    words = {
        "hit": "reused from the cache",
        "miss": "recomputed from scratch",
        "auto": (
            "reused wherever the old run actually holds the same job, and "
            "recomputed otherwise"
        ),
        "ignore": "not checked",
    }
    lines.append(f"By default every CNS output should be {words.get(default, default)}.")
    for module, verdict in (expect.get("modules") or {}).items():
        lines.append(f"  ... except {module} output, which should be {words[verdict]}.")
    for rule in expect.get("patterns") or []:
        lines.append(
            f"  ... except files matching {rule['match']}, "
            f"which should be {words[rule['expect']]}."
        )
    for path, source in (expect.get("paths") or {}).items():
        if source is None:
            lines.append(f"  ... except {path}, which should be recomputed.")
        else:
            lines.append(f"  ... except {path}, which should be reused from {source}.")
    return "\n".join(lines)


def _describe_perturbation(case: Case) -> list[str]:
    """One readable sentence per perturbation operation."""
    sentences = []
    for operation in case.perturb:
        op = operation.get("op")
        if op == "top":
            for key, value in (operation.get("set") or {}).items():
                sentences.append(f"set {key} = {value!r} at the top of the config")
        elif op == "module":
            module = operation["module"]
            for key, value in (operation.get("set") or {}).items():
                sentences.append(f"set {key} = {value!r} in [{module}]")
            for key in operation.get("unset") or []:
                sentences.append(f"remove {key} from [{module}]")
        elif op == "insert":
            sentences.append(f"add a [{operation['module']}] step")
        elif op == "remove":
            sentences.append(f"delete the [{operation['module']}] step")
        elif op == "duplicate":
            sentences.append(f"repeat the [{operation['module']}] step")
        elif op == "swap_steps":
            sentences.append(f"swap the [{operation['a']}] and [{operation['b']}] steps")
        elif op == "molecules":
            if operation.get("swap"):
                sentences.append("swap two entries in the molecules list")
            elif operation.get("reverse"):
                sentences.append("reverse the molecules list")
        elif op == "edit_input":
            how = {
                "coordinate": "move one atom by 1 Angstrom in",
                "remark": "add a REMARK line to",
                "whitespace": "change the whitespace and line endings of",
                "segid": "rename a chain throughout",
                "restraint": "change one distance in",
                "gzip": "gzip",
                "models": "change the ensemble members of",
            }.get(operation["how"], operation["how"])
            sentences.append(f"{how} {operation['file']}")
        elif op == "copy_input":
            sentences.append(
                f"copy {operation['file']} to {operation['as']} and use the copy "
                "(same content, different name)"
            )
        elif op == "move_input":
            sentences.append(f"move {operation['file']} somewhere else")
        elif op == "env":
            for key, value in (operation.get("set") or {}).items():
                sentences.append(f"set the environment variable {key}={value}")
        elif op == "cwd":
            sentences.append("run haddock3 from a different working directory")
        elif op == "install":
            edits = operation.get("edit") or {}
            if edits:
                for relative in edits:
                    sentences.append(f"edit the installed HADDOCK3 file {relative}")
            else:
                sentences.append("run against a second copy of HADDOCK3 at another path")
        elif op == "cns_wrapper":
            sentences.append("point cns_exec at a shell wrapper around the real CNS")
        elif op == "cosmetic":
            sentences.append("reorder the config parameters and add comments")
        elif op == "explicit_default":
            sentences.append(
                f"write {operation['param']} out explicitly at its default value "
                f"in [{operation['module']}]"
            )
        elif op == "restraints_archive":
            sentences.append(
                f"give every model its own restraint table "
                f"({operation.get('members', 'identical')} ones)"
            )
        elif op == "note":
            sentences.append(operation.get("text", ""))
    return [s for s in sentences if s]


def explain(case: Case) -> str:
    out = [f"{case.case}", "=" * len(case.case), ""]
    if case.title:
        out += [f"  {case.title}", ""]
    out += [f"  Taxonomy: {case.taxonomy or case.axis}", ""]
    if case.why:
        out += ["Why this matters", "----------------", _wrap(case.why), ""]
    out += [
        "The two runs",
        "------------",
        f"  A (the old run):  the '{case.base}' workflow in the corpus",
        f"  cache source(s):  {', '.join(case.sources) or '(none)'}",
        "",
        "  B is A with these changes:",
    ]
    changes = _describe_perturbation(case)
    out += [f"    - {sentence}" for sentence in changes] or ["    - nothing at all"]
    out += ["", "What should happen", "------------------", _plain_expectation(case), ""]
    if case.mode != "complete":
        out += [
            f"  This case runs in '{case.mode}' mode rather than running to "
            "completion; see the suite README for what that means.",
            "",
        ]
    if case.degraded:
        out += ["Note", "----", _wrap(case.degraded), ""]
    return "\n".join(out)


RUN_SH = """#!/bin/bash
# Everything this test does, as two ordinary HADDOCK3 commands.
#
# Run it, then edit B.cfg and run it again -- the last command prints which
# files were reused and which were recomputed, so you can see immediately what
# your edit cost.
set -e
cd "$(dirname "$0")"
{env}
echo "== Run A (the 'before' run; this is the one that costs CNS time) =="
rm -rf runA
haddock3 A.cfg

echo
echo "== Run B (the 'after' run; --cache lets it reuse A's results) =="
rm -rf runB
haddock3 B.cfg --cache runA

echo
echo "== What was reused =="
python {checker} runA runB
"""


def export(case: Case, destination: Path, corpus: Manifest, corpus_dir: Path) -> None:
    """Write a self-contained folder that needs nothing from this suite."""
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    work = destination / ".build"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()

    prepared = prepare(case, corpus, corpus_dir, work)

    # The 'before' run: the base workflow itself, run locally rather than
    # taken from the corpus, so the exported folder stands alone.
    # A's inputs are the corpus's; B's are either the same or a private copy
    # the perturbation edited. Export both, so an edited input is visible as a
    # file you can diff rather than as something the harness did.
    shutil.copytree(
        corpus_dir / "systems" / _system(case), destination / "data", dirs_exist_ok=True
    )
    b_data = "data"
    if prepared.data_dir != corpus_dir / "systems" / _system(case):
        shutil.copytree(prepared.data_dir, destination / "data-B", dirs_exist_ok=True)
        b_data = "data-B"

    base_config = (corpus_dir / "configs" / f"{case.base}.cfg").read_text(
        encoding="utf-8"
    )
    base_config = base_config.replace(
        str(corpus_dir / "systems" / _system(case)), "data"
    )
    base_config = _set_run_dir(base_config, "runA")
    (destination / "A.cfg").write_text(base_config, encoding="utf-8")

    b_config = prepared.config_path.read_text(encoding="utf-8")
    b_config = b_config.replace(str(prepared.data_dir), b_data)
    b_config = b_config.replace(str(corpus_dir / "systems" / _system(case)), "data")
    b_config = b_config.replace(str(work), ".")
    b_config = _set_run_dir(b_config, "runB")
    (destination / "B.cfg").write_text(b_config, encoding="utf-8")

    env_lines = []
    for key, value in prepared.env.items():
        env_lines.append(f"export {key}={value}")
    if prepared.install is not None:
        shutil.copytree(prepared.install, destination / "install", dirs_exist_ok=True)
        env_lines.append('export PYTHONPATH="$PWD/install:$PYTHONPATH"')
    shutil.copy2(CHECK_REUSE, destination / "check_reuse.py")
    (destination / "run.sh").write_text(
        RUN_SH.format(env="\n".join(env_lines), checker="check_reuse.py"),
        encoding="utf-8",
    )
    (destination / "run.sh").chmod(0o755)

    readme = [
        f"# {case.title or case.case}",
        "",
        explain(case),
        "",
        "## How to run this",
        "",
        "```bash",
        "./run.sh",
        "```",
        "",
        "`run.sh` does three things: it runs `A.cfg` (the 'before' workflow),",
        "then runs `B.cfg` with `--cache runA` (the 'after' workflow, allowed to",
        "reuse A's results), then prints which files were reused. Everything in",
        "this folder is self-contained: `check_reuse.py` imports nothing but the",
        "Python standard library.",
        "",
        "## Checking it by hand, with no scripts at all",
        "",
        "A reused result is a *hardlink*: the old file and the new one are the",
        "same file on disk. `ls -i` prints the inode number, and two names with",
        "the same inode number are the same file:",
        "",
        "```bash",
        "ls -i runA/*/rigidbody_1.pdb runB/*/rigidbody_1.pdb",
        "```",
        "",
        "Same number: reused. Different numbers: recomputed. That is the whole",
        "measurement, and `check_reuse.py` is only a tidier way of doing it for",
        "every file at once.",
        "",
        "## Things worth trying",
        "",
        "* Edit `B.cfg` and rerun. Any change that alters what CNS computes",
        "  should turn some files from REUSED into RECOMPUTED; any change that",
        "  only alters names, paths or bookkeeping should not.",
        "* Rename `runB` in `B.cfg` to something else. Nothing should change:",
        "  where a run lives is not part of what its results are.",
        "* Delete a file from `runA` and rerun B. That one file should be",
        "  recomputed and the rest still reused.",
        "",
    ]
    if case.sources != [case.base]:
        readme += [
            "## Note on the cache source",
            "",
            "This case uses a specially prepared cache source "
            f"(`{', '.join(case.sources)}`) rather than a plain previous run.",
            "`run.sh` above uses a plain `runA` instead, so it demonstrates the",
            "shape of the case but not the specially prepared part. The full",
            "version lives in the corpus; see the suite README.",
            "",
        ]
    (destination / "README.md").write_text("\n".join(readme), encoding="utf-8")
    shutil.rmtree(work, ignore_errors=True)


def _system(case: Case) -> str:
    from cachesuite.systems import BASE_RUNS

    return next(spec.system for spec in BASE_RUNS if spec.name == case.base)


def _set_run_dir(text: str, name: str) -> str:
    out = []
    for line in text.splitlines():
        if line.strip().startswith("run_dir"):
            out.append(f'run_dir = "{name}"')
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def run_once(case: Case, corpus: Manifest, corpus_dir: Path, workdir: Path) -> int:
    """Run the case once and print the reuse table.  No gates, no verdict."""
    from cachesuite.cases import execute

    workdir.mkdir(parents=True, exist_ok=True)
    prepared = prepare(case, corpus, corpus_dir, workdir)
    print(f"running: {' '.join(prepared.config_path.parts[-1:])} "
          f"--cache {' --cache '.join(str(s) for s in prepared.sources.values())}")
    result = execute(prepared, timeout=3600)
    print(f"exit code {result.returncode}, {result.duration:.1f} s")
    if not result.run_dir.is_dir():
        print(result.tail(30))
        return 1
    print()
    sources = list(prepared.sources.values())
    artifacts = cacheable_artifacts(result.run_dir)
    width = max((len(a.relative) for a in artifacts), default=10)
    for artifact in artifacts:
        origin = link_source(result.run_dir / artifact.relative, sources)
        if origin is None:
            print(f"{artifact.relative:<{width}}  RECOMPUTED")
        else:
            print(f"{artifact.relative:<{width}}  REUSED       <- {origin.name}")
    print()
    print("Compare that against 'What should happen' above.")
    return 0 if result.ok else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("case", nargs="?", help="the case name")
    parser.add_argument("--list", action="store_true", help="list the cases")
    parser.add_argument(
        "--common",
        action="store_true",
        help="with --list, show only the everyday cases",
    )
    parser.add_argument("--run", action="store_true", help="run the case once")
    parser.add_argument("--export", type=Path, help="write a standalone folder")
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("./try-case-work"),
        help="where --run puts its files (default: ./try-case-work)",
    )
    arguments = parser.parse_args(argv)

    cases = {case.case: case for case in load_cases()}

    if arguments.list or not arguments.case:
        for case in cases.values():
            if arguments.common and not case.common:
                continue
            if case.skip:
                continue
            marker = "*" if case.common else " "
            print(f"{marker} {case.case}")
            print(f"    {case.title}")
        print()
        print("* = an everyday case, and a good place to start.")
        print("Use: try_case.py <name> [--run | --export DIR]")
        return 0

    case = cases.get(arguments.case)
    if case is None:
        matches = [name for name in cases if arguments.case in name]
        parser.error(
            f"no case named {arguments.case!r}"
            + (f"; did you mean one of {matches}?" if matches else "")
        )

    print(explain(case))
    if not (arguments.run or arguments.export):
        return 0
    if case.skip:
        print("This case is deliberately not run; there is nothing to execute.")
        return 0

    corpus_dir = corpus_root()
    try:
        corpus = Manifest.read(corpus_dir)
    except CorpusNotBuilt as error:
        print(f"\n{error}")
        return 1

    if arguments.export:
        export(case, arguments.export, corpus, corpus_dir)
        print(f"\nWrote a standalone copy to {arguments.export.resolve()}")
        print("Run it with:  cd", arguments.export, "&& ./run.sh")
        return 0
    return run_once(case, corpus, corpus_dir, arguments.workdir / case.case)


if __name__ == "__main__":
    raise SystemExit(main())
