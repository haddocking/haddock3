"""
Run the example test cases and compare them with a reference from `main`.

Every `examples/*/*-test.cfg` is run, except MPI/batch examples and follow-up
configurations that need an existing run (see `skip_reason`). For each run the
script checks that it completed without errors and, for every caprieval step,
extracts the CAPRI statistics also reported by `extract-capri-stats-dockq.sh`
(and `extract-capri-stats-ligand.sh` for protein-ligand examples).

The reference is built from `origin/main` in a separate git worktree, so the
current checkout is never touched. All results are written to
`examples/regression/`.

WARNING: existing run directories with the same name as those defined in the
test configuration files are deleted before each run.

USAGE:

    $ python run_examples_regression.py init      # build the main reference
    $ python run_examples_regression.py run       # run current code + compare
    $ python run_examples_regression.py compare   # only compare saved results
    $ python run_examples_regression.py --list    # show selected/skipped tests
    $ python run_examples_regression.py run --only scoring/ --only protein-DNA
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

try:
    from haddock.gear.config import load as read_config
except Exception:  # noqa: BLE001
    print(
        "Haddock3 could not be imported. "
        "Please activate the haddock3 python environment.",
        file=sys.stderr,
    )
    sys.exit(1)


EXAMPLES = Path(__file__).resolve().parent
REPO = EXAMPLES.parent
OUTDIR = Path(EXAMPLES, "regression")
REFERENCE_JSON = Path(OUTDIR, "reference.json")
CURRENT_JSON = Path(OUTDIR, "current.json")
REPORT_MD = Path(OUTDIR, "report.md")
DEFAULT_WORKTREE = Path(REPO.parent, "haddock3-main-ref")

# configuration files that need an existing run (--extend-run / --restart)
FOLLOW_UPS = {
    "docking-extend-run-exit-test.cfg": "needs an existing run (--extend-run)",
    "docking-extend-run-test.cfg": "needs an existing run (--extend-run)",
    "docking-restart-exit-test.cfg": "needs an existing run (--restart)",
}

# CAPRI quality thresholds: (column, acceptable, medium, high, higher_is_better)
DOCKQ_CRITERIA = ("dockq", 0.23, 0.49, 0.80, True)
LIGAND_CRITERIA = ("ilrmsd", 2.5, 1.0, 0.5, False)
LIGAND_FOLDERS = ("docking-protein-ligand",)

FINISHED_LINE = "cli INFO] Finished at"
ERROR_LINE = re.compile(r"^\[[^\]]* (ERROR|CRITICAL)\]")
STEP_FOLDER = re.compile(r"^(\d+)_(\w+)$")
IGNORED_SUFFIXES = (".out", ".err", ".inp", ".seed")


def printf(*args):
    """Print and flush."""
    print(*args, flush=True)


# ---------------------------------------------------------------------------
# Test selection
# ---------------------------------------------------------------------------


def discover_tests(examples_dir):
    """Return a list of (test id, cfg path, skip reason or None)."""
    tests = []
    for cfg in sorted(examples_dir.glob("*/*-test.cfg")):
        test_id = f"{cfg.parent.name}/{cfg.stem}"
        tests.append((test_id, cfg, skip_reason(cfg)))
    return tests


def skip_reason(cfg):
    """Return why `cfg` cannot be run here, or None."""
    if "mpi" in cfg.name:
        return "MPI example"
    if cfg.name in FOLLOW_UPS:
        return FOLLOW_UPS[cfg.name]
    try:
        params = read_config(cfg)["final_cfg"]
    except Exception as err:  # noqa: BLE001
        return f"configuration could not be read ({err})"
    if params.get("mode", "local") != "local":
        return f"mode = {params['mode']}"
    return None


def select(tests, patterns):
    """Keep tests whose id contains any of `patterns`."""
    if not patterns:
        return tests
    return [t for t in tests if any(p in t[0] for p in patterns)]


def expected_steps(cfg):
    """Return the module names expected as step folders, in order."""
    params = read_config(cfg)["final_cfg"]
    modules = []
    for key, value in params.items():
        if isinstance(value, dict) and "." in key:
            name = key.rsplit(".", 1)[0]
            if name == "exit":
                break
            modules.append(name)
    return modules


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------


def run_test(test_id, cfg, env, log_dir, timeout):
    """Run one example and return its result dictionary."""
    run_dir = Path(cfg.parent, read_config(cfg)["final_cfg"]["run_dir"])
    shutil.rmtree(run_dir, ignore_errors=True)

    log_file = Path(log_dir, test_id.replace("/", "__") + ".log")
    log_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        "-c",
        "from haddock.clis.cli import maincli; maincli()",
        cfg.name,
    ]
    start = time.time()
    with open(log_file, "w") as fout:
        try:
            proc = subprocess.run(
                cmd,
                cwd=cfg.parent,
                env=env,
                stdout=fout,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
            returncode = proc.returncode
        except subprocess.TimeoutExpired:
            returncode = "timeout"
    runtime = time.time() - start

    result = {
        "returncode": returncode,
        "runtime": round(runtime, 1),
        "log_tail": log_file.read_text(errors="replace").splitlines()[-15:],
    }
    result.update(inspect_run(run_dir, expected_steps(cfg), cfg.parent.name))
    result["status"] = "OK" if not result["problems"] else "FAILED"
    if returncode != 0:
        result["problems"].insert(0, f"exit code {returncode}")
        result["status"] = "FAILED"
    return result


def inspect_run(run_dir, modules, folder):
    """Check a finished run directory and extract its statistics."""
    problems = []
    log = Path(run_dir, "log")
    errors = []
    if not log.exists():
        problems.append("no run log")
    else:
        lines = log.read_text(errors="replace").splitlines()
        errors = [ln for ln in lines if ERROR_LINE.match(ln)]
        if not any(FINISHED_LINE in ln for ln in lines):
            problems.append("run did not finish ('Finished at' missing in log)")
        if errors:
            problems.append(f"{len(errors)} ERROR/CRITICAL line(s) in log")

    steps = {}
    if run_dir.is_dir():
        for path in run_dir.iterdir():
            match = STEP_FOLDER.match(path.name)
            if path.is_dir() and match:
                steps[int(match.group(1))] = path

    for index, module in enumerate(modules):
        path = steps.get(index)
        if path is None or STEP_FOLDER.match(path.name).group(2) != module:
            problems.append(f"step {index}_{module} missing")

    criteria = [DOCKQ_CRITERIA]
    if folder.startswith(LIGAND_FOLDERS):
        criteria.append(LIGAND_CRITERIA)

    step_info = {}
    for index in sorted(steps):
        path = steps[index]
        info = fingerprint(path)
        capri = Path(path, "capri_ss.tsv")
        if capri.exists():
            info["capri"] = {c[0]: capri_stats(capri, c) for c in criteria}
        # key without zero padding, so 01_x and 1_x match between runs
        step_info[f"{index}_{STEP_FOLDER.match(path.name).group(2)}"] = info

    return {"problems": problems, "errors": errors[:20], "steps": step_info}


def fingerprint(step_dir):
    """Return the number of output models and the output file names."""
    io_json = Path(step_dir, "io.json")
    nmodels = None
    if io_json.exists():
        try:
            nmodels = len(json.loads(io_json.read_text())["output"])
        except Exception:  # noqa: BLE001
            nmodels = None
    files = sorted(
        {
            p.name.removesuffix(".gz")
            for p in step_dir.iterdir()
            if not p.name.removesuffix(".gz").endswith(IGNORED_SUFFIXES)
        }
    )
    return {"io_json": io_json.exists(), "nmodels": nmodels, "files": files}


def _float(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


def capri_stats(capri_file, criteria):
    """Extract CAPRI statistics from a capri_ss.tsv file."""
    column, accept, medium, high, higher_better = criteria
    with open(capri_file) as fin:
        rows = list(csv.DictReader(fin, delimiter="\t"))

    def passes(row, cutoff):
        value = _float(row.get(column))
        if value is None:
            return False
        return value >= cutoff if higher_better else value <= cutoff

    def summary(row):
        if row is None:
            return None
        return {
            "rank": row.get("caprieval_rank"),
            "irmsd": _float(row.get("irmsd")),
            "ilrmsd": _float(row.get("ilrmsd")),
            "fnat": _float(row.get("fnat")),
            "dockq": _float(row.get("dockq")),
        }

    # best model: lowest i-RMSD (DockQ criteria) or il-RMSD (ligand criteria)
    best_col = "irmsd" if column == "dockq" else "ilrmsd"
    ranked = [r for r in rows if _float(r.get(best_col)) is not None]
    best = min(ranked, key=lambda r: _float(r[best_col])) if ranked else None

    return {
        "nmodels": len(rows),
        "acceptable": sum(passes(r, accept) for r in rows),
        "medium": sum(passes(r, medium) for r in rows),
        "high": sum(passes(r, high) for r in rows),
        "first_acceptable": summary(next((r for r in rows if passes(r, accept)), None)),
        "first_medium": summary(next((r for r in rows if passes(r, medium)), None)),
        "best": summary(best),
    }


def run_all(tests, env, log_dir, timeout, break_on_errors):
    """Run the selected tests and return the results."""
    results = {}
    runnable = [t for t in tests if t[2] is None]
    for i, (test_id, cfg, _) in enumerate(runnable, start=1):
        printf(f"[{i}/{len(runnable)}] {test_id} ...")
        try:
            result = run_test(test_id, cfg, env, log_dir, timeout)
        except Exception as err:  # noqa: BLE001
            # never let one broken example stop a long series of runs
            result = {
                "status": "FAILED",
                "returncode": None,
                "runtime": 0.0,
                "problems": [f"could not be run: {err!r}"],
                "errors": [],
                "log_tail": [],
                "steps": {},
            }
        results[test_id] = result
        printf(f"    {result['status']} ({result['runtime']:.0f} s)")
        for problem in result["problems"]:
            printf(f"      - {problem}")
        if break_on_errors and result["status"] != "OK":
            break
    return results


def save(path, label, commit, tests, results, previous=None):
    """Save the results, merging them with a previous file for --only runs."""
    data = previous or {"tests": {}, "skipped": {}}
    data.update(
        {
            "label": label,
            "commit": commit,
            "date": datetime.now().isoformat(timespec="seconds"),
        }
    )
    data["tests"].update(results)
    data["skipped"].update({t[0]: t[2] for t in tests if t[2] is not None})
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1))
    printf(f"Results saved to {path}")


def load_json(path):
    """Load a results file if it exists."""
    return json.loads(path.read_text()) if path.exists() else None


def git(*args, cwd=REPO):
    """Run a git command and return its stdout."""
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


def prepare_worktree(worktree, ref):
    """Create or update a detached worktree at `ref`."""
    git("fetch", "origin")
    if not Path(worktree, ".git").exists():
        git("worktree", "add", "--detach", str(worktree), ref)
    else:
        git("checkout", "--detach", ref, cwd=worktree)
    # compiled helpers are not tracked by git
    shutil.copytree(
        Path(REPO, "src/haddock/bin"),
        Path(worktree, "src/haddock/bin"),
        dirs_exist_ok=True,
    )
    return git("rev-parse", "--short", "HEAD", cwd=worktree)


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------


def _fmt(value, digits=3):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def _arrow(old, new, higher_better):
    if old is None or new is None or old == new:
        return ""
    better = new > old if higher_better else new < old
    return " ▲" if better else " ▼"


def compare_capri(ref, new, tol_dockq, tol_irmsd):
    """Return a list of differences between two capri_stats dictionaries."""
    diffs = []
    for key in ("nmodels", "acceptable", "medium", "high"):
        if ref[key] != new[key]:
            arrow = "" if key == "nmodels" else _arrow(ref[key], new[key], True)
            diffs.append(f"{key}: {ref[key]} → {new[key]}{arrow}")

    def rank(stats, key):
        return stats[key]["rank"] if stats[key] else None

    for key in ("first_acceptable", "first_medium"):
        old_rank, new_rank = rank(ref, key), rank(new, key)
        if old_rank != new_rank:
            arrow = ""
            if old_rank is not None and new_rank is not None:
                arrow = _arrow(int(old_rank), int(new_rank), False)
            diffs.append(
                f"{key.replace('_', ' ')} rank: {_fmt(old_rank)} → "
                f"{_fmt(new_rank)}{arrow}"
            )

    old_best, new_best = ref["best"] or {}, new["best"] or {}
    if old_best.get("rank") != new_best.get("rank"):
        diffs.append(
            f"best model rank: {_fmt(old_best.get('rank'))} → "
            f"{_fmt(new_best.get('rank'))}"
            f"{_arrow(_int(old_best.get('rank')), _int(new_best.get('rank')), False)}"
        )
    for key, tol, higher in (
        ("dockq", tol_dockq, True),
        ("irmsd", tol_irmsd, False),
        ("ilrmsd", tol_irmsd, False),
    ):
        old, cur = old_best.get(key), new_best.get(key)
        if old is None or cur is None:
            continue
        if abs(cur - old) > tol:
            diffs.append(
                f"best {key}: {_fmt(old)} → {_fmt(cur)}{_arrow(old, cur, higher)}"
            )
    return diffs


def _int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _rank_change(old, new):
    """Return -1 (worse), 0 (same) or 1 (better) for a rank change.

    A missing rank (no such model) is worse than any rank.
    """
    old, new = _int(old), _int(new)
    if old == new:
        return 0
    if old is None:
        return 1
    if new is None:
        return -1
    return 1 if new < old else -1


def _value_change(old, new, tol, higher_better):
    """Return -1, 0 or 1 for a value change beyond `tol`."""
    if old is None or new is None or abs(new - old) <= tol:
        return 0
    return 1 if (new > old) == higher_better else -1


def _cell(old, new, change, digits=3):
    """Format 'old → new ▲/▼', or just the value when unchanged."""
    if _fmt(old, digits) == _fmt(new, digits):
        return _fmt(new, digits)
    arrow = {1: " ▲", -1: " ▼"}.get(change, "")
    return f"{_fmt(old, digits)} → {_fmt(new, digits)}{arrow}"


def ranking_row(test_id, step, crit, ref, new, tol_dockq, tol_irmsd):
    """Analyse the first acceptable and the best model of one caprieval step.

    Returns the markdown table row and the changes as a dictionary:
    `first` (rank of the first acceptable model), `best_rank` and
    `best_quality` (DockQ, or il-RMSD for the ligand criteria), each -1, 0, 1.
    """
    old_first = (ref["first_acceptable"] or {}).get("rank")
    new_first = (new["first_acceptable"] or {}).get("rank")
    old_best, new_best = ref["best"] or {}, new["best"] or {}

    if crit == "dockq":
        qkey, tol, higher = "dockq", tol_dockq, True
        rmsd_key = "irmsd"
    else:
        qkey, tol, higher = "ilrmsd", tol_irmsd, False
        rmsd_key = "ilrmsd"

    changes = {
        "first": _rank_change(old_first, new_first),
        "best_rank": _rank_change(old_best.get("rank"), new_best.get("rank")),
        "best_quality": _value_change(
            old_best.get(qkey), new_best.get(qkey), tol, higher
        ),
    }
    rmsd_change = _value_change(
        old_best.get(rmsd_key), new_best.get(rmsd_key), tol_irmsd, False
    )
    dockq_change = _value_change(
        old_best.get("dockq"), new_best.get("dockq"), tol_dockq, True
    )

    values = set(changes.values())
    if values == {0}:
        verdict = "same"
    elif -1 not in values:
        verdict = "**better**"
    elif 1 not in values:
        verdict = "**worse**"
    else:
        verdict = "mixed"

    label = "DockQ" if crit == "dockq" else "ligand"
    row = (
        f"| {test_id} | {step} | {label} "
        f"| {_cell(old_first, new_first, changes['first'])} "
        f"| {_cell(old_best.get('rank'), new_best.get('rank'), changes['best_rank'])} "
        f"| {_cell(old_best.get(rmsd_key), new_best.get(rmsd_key), rmsd_change)} "
        f"| {_cell(old_best.get('fnat'), new_best.get('fnat'), 0)} "
        f"| {_cell(old_best.get('dockq'), new_best.get('dockq'), dockq_change)} "
        f"| {verdict} |"
    )
    return row, changes


def ranking_section(rows, tallies):
    """Return the report lines of the first-acceptable/best-model analysis."""
    if not rows:
        return ["No caprieval steps to compare."]

    def tally(key):
        counts = [tallies[key].count(v) for v in (1, 0, -1)]
        return f"{counts[0]} better, {counts[1]} unchanged, {counts[2]} worse"

    header = [
        (
            "| Test | Step | Criteria | First acceptable rank | Best model rank "
            "| Best i-RMSD / il-RMSD | Best Fnat | Best DockQ | Verdict |"
        ),
        "|---|---|---|---|---|---|---|---|---|",
    ]
    changed = [row for row in rows if not row.endswith("| same |")]
    return [
        f"Over {len(rows)} caprieval step(s):",
        "",
        f"- rank of the first acceptable or better model: {tally('first')}",
        f"- rank of the best model: {tally('best_rank')}",
        (
            "- quality of the best model (DockQ, or il-RMSD for ligands): "
            f"{tally('best_quality')}"
        ),
        "",
        (
            "The best model is the one with the lowest i-RMSD (il-RMSD for the "
            "ligand criteria). Values are `main → now`; a single value means "
            "unchanged."
        ),
        "",
        f"### Steps that changed ({len(changed)})",
        "",
        *((header + changed) if changed else ["None."]),
        "",
        "<details>",
        f"<summary>All caprieval steps ({len(rows)})</summary>",
        "",
        *header,
        *rows,
        "",
        "</details>",
    ]


def compare_steps(ref_steps, new_steps):
    """Compare step fingerprints; return (output problems, capri pairs)."""
    problems = []
    capri_pairs = []
    for step, ref_info in ref_steps.items():
        new_info = new_steps.get(step)
        if new_info is None:
            problems.append(f"{step}: step folder missing")
            continue
        if ref_info["io_json"] and not new_info["io_json"]:
            problems.append(f"{step}: io.json missing")
        if (
            ref_info["nmodels"] is not None
            and new_info["nmodels"] is not None
            and new_info["nmodels"] < ref_info["nmodels"]
        ):
            problems.append(
                f"{step}: {new_info['nmodels']} output models "
                f"(reference {ref_info['nmodels']})"
            )
        missing = sorted(set(ref_info["files"]) - set(new_info["files"]))
        if missing:
            shown = ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
            problems.append(f"{step}: {len(missing)} file(s) missing: {shown}")
        if "capri" in ref_info and "capri" in new_info:
            capri_pairs.append((step, ref_info["capri"], new_info["capri"]))
    return problems, capri_pairs


def compare(ref, cur, tol_dockq, tol_irmsd, tol_runtime):
    """Compare two result files and return (report lines, summary rows)."""
    status_lines, quality_lines, runtime_lines = [], [], []
    failed_both, unchanged = [], []
    ranking_rows = []
    ranking_tallies = {"first": [], "best_rank": [], "best_quality": []}
    summary = []

    for test_id in sorted(set(ref["tests"]) | set(cur["tests"])):
        r, c = ref["tests"].get(test_id), cur["tests"].get(test_id)
        if r is None or c is None:
            where = "reference" if r is None else "current run"
            status_lines.append(f"- **{test_id}**: missing in the {where}")
            summary.append((test_id, "MISSING"))
            continue

        if r["status"] != "OK" and c["status"] != "OK":
            failed_both.append(f"- {test_id}: {'; '.join(c['problems'])}")
            summary.append((test_id, "FAILED (also on main)"))
            continue

        if c["status"] != "OK":
            status_lines.append(f"- **{test_id}**: OK on main, **FAILED** now")
            status_lines += [f"  - {p}" for p in c["problems"]]
            status_lines += [f"  - `{e}`" for e in c["errors"][:5]]
            if c["log_tail"]:
                status_lines.append("  - last log lines:")
                status_lines.append("    ```")
                status_lines += [f"    {ln}" for ln in c["log_tail"]]
                status_lines.append("    ```")
            summary.append((test_id, "FAILED"))
            continue

        changed = False
        if r["status"] != "OK":
            status_lines.append(f"- **{test_id}**: FAILED on main, OK now")
            summary.append((test_id, "FIXED"))
            continue

        output_problems, capri_pairs = compare_steps(r["steps"], c["steps"])
        if output_problems:
            changed = True
            status_lines.append(f"- **{test_id}**: completed, but outputs differ")
            status_lines += [f"  - {p}" for p in output_problems]

        for step, ref_capri, new_capri in capri_pairs:
            for crit, ref_stats in ref_capri.items():
                if crit not in new_capri:
                    continue
                row, changes = ranking_row(
                    test_id,
                    step,
                    crit,
                    ref_stats,
                    new_capri[crit],
                    tol_dockq,
                    tol_irmsd,
                )
                ranking_rows.append(row)
                for key, value in changes.items():
                    ranking_tallies[key].append(value)
                diffs = compare_capri(ref_stats, new_capri[crit], tol_dockq, tol_irmsd)
                if diffs:
                    changed = True
                    label = "DockQ" if crit == "dockq" else "ligand il-RMSD"
                    quality_lines.append(
                        f"- **{test_id}** `{step}` ({label} criteria): "
                        + "; ".join(diffs)
                    )

        if r["runtime"] > 0:
            ratio = c["runtime"] / r["runtime"]
            if abs(ratio - 1) > tol_runtime:
                runtime_lines.append(
                    f"- {test_id}: {r['runtime']:.0f} s → {c['runtime']:.0f} s "
                    f"({(ratio - 1) * 100:+.0f}%)"
                )

        if changed:
            summary.append((test_id, "CHANGED"))
        else:
            unchanged.append(test_id)
            summary.append((test_id, "OK"))

    report = [
        "# Example regression report",
        "",
        f"- Reference: {ref.get('label')} `{ref.get('commit')}` ({ref.get('date')})",
        f"- Current: {cur.get('label')} `{cur.get('commit')}` ({cur.get('date')})",
        (
            f"- Tolerances: DockQ {tol_dockq}, i-RMSD/il-RMSD {tol_irmsd} Å, "
            f"runtime {tol_runtime * 100:.0f}%"
        ),
        "- ▲ better, ▼ worse (higher DockQ and lower (i)RMSD / rank are better)",
        "",
        "## 1. Run status and outputs",
        *(status_lines or ["No status or output changes."]),
        "",
        "## 2. CAPRI quality changes",
        *(quality_lines or ["No changes in CAPRI statistics."]),
        "",
        "## 3. First acceptable and best models",
        *ranking_section(ranking_rows, ranking_tallies),
        "",
        "## 4. Runtime changes",
        *(runtime_lines or ["No significant runtime changes."]),
        "",
        "## 5. Failing on both main and the current code",
        *(failed_both or ["None."]),
        "",
        "## 6. Unchanged",
        *([f"- {t}" for t in unchanged] or ["None."]),
        "",
        "## 7. Skipped",
        *[
            f"- {t}: {why}"
            for t, why in sorted({**ref["skipped"], **cur["skipped"]}.items())
        ],
        "",
    ]
    return report, summary


def write_report(args):
    """Compare reference and current results and write the report."""
    ref, cur = load_json(REFERENCE_JSON), load_json(CURRENT_JSON)
    if ref is None or cur is None:
        missing = REFERENCE_JSON if ref is None else CURRENT_JSON
        printf(f"{missing} not found: run `init` and `run` first.")
        sys.exit(1)
    report, summary = compare(
        ref, cur, args.tol_dockq, args.tol_irmsd, args.tol_runtime
    )
    REPORT_MD.write_text(os.linesep.join(report))
    printf("")
    width = max(len(t) for t, _ in summary) if summary else 10
    for test_id, state in summary:
        printf(f"{test_id:<{width}}  {state}")
    printf(f"\nReport written to {REPORT_MD}")


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def load_args():
    """Parse the command line."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "command",
        nargs="?",
        choices=("init", "run", "compare"),
        help="init: run the reference from main; run: run the current code "
        "and compare; compare: only compare saved results.",
    )
    ap.add_argument(
        "--list",
        action="store_true",
        help="List the selected and skipped tests, then exit.",
    )
    ap.add_argument(
        "--only",
        action="append",
        default=[],
        metavar="PATTERN",
        help="Only run tests whose id (folder/cfg-name) contains PATTERN. "
        "Can be repeated.",
    )
    ap.add_argument(
        "--ref-worktree",
        type=Path,
        default=DEFAULT_WORKTREE,
        help=f"Worktree used for the reference (default: {DEFAULT_WORKTREE}).",
    )
    ap.add_argument(
        "--ref-branch",
        default="origin/main",
        help="Git reference used to build the reference (default: origin/main).",
    )
    ap.add_argument(
        "--timeout",
        type=float,
        default=7200,
        help="Maximum time per test in seconds (default: 7200).",
    )
    ap.add_argument(
        "-b",
        "--break-on-errors",
        action="store_true",
        help="Stop as soon as a test fails.",
    )
    ap.add_argument("--tol-dockq", type=float, default=0.05)
    ap.add_argument("--tol-irmsd", type=float, default=0.5)
    ap.add_argument(
        "--tol-runtime",
        type=float,
        default=0.5,
        help="Relative runtime change reported (default: 0.5 = 50%%).",
    )
    args = ap.parse_args()
    if not args.list and args.command is None:
        ap.error("a command (init, run or compare) or --list is required")
    return args


def list_tests(tests):
    """Print the selected and skipped tests."""
    for test_id, cfg, reason in tests:
        if reason:
            printf(f"SKIP  {test_id}  ({reason})")
            continue
        run_dir = Path(cfg.parent, read_config(cfg)["final_cfg"]["run_dir"])
        note = "  [existing run dir will be deleted]" if run_dir.exists() else ""
        printf(f"RUN   {test_id}{note}")
    nrun = sum(1 for t in tests if t[2] is None)
    printf(f"\n{nrun} tests selected, {len(tests) - nrun} skipped.")


def main():
    """Run the command."""
    args = load_args()
    tests = select(discover_tests(EXAMPLES), args.only)

    if args.list:
        list_tests(tests)
        return

    if args.command == "compare":
        write_report(args)
        return

    env = dict(os.environ)
    if args.command == "init":
        worktree = args.ref_worktree.resolve()
        commit = prepare_worktree(worktree, args.ref_branch)
        printf(f"Reference worktree: {worktree} at {args.ref_branch} ({commit})")
        env["PYTHONPATH"] = os.pathsep.join(
            filter(None, [str(Path(worktree, "src")), env.get("PYTHONPATH")])
        )
        ref_tests = select(discover_tests(Path(worktree, "examples")), args.only)
        results = run_all(
            ref_tests,
            env,
            Path(OUTDIR, "logs", "reference"),
            args.timeout,
            args.break_on_errors,
        )
        previous = load_json(REFERENCE_JSON) if args.only else None
        save(REFERENCE_JSON, args.ref_branch, commit, ref_tests, results, previous)
        return

    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    commit = git("rev-parse", "--short", "HEAD")
    if git("status", "--porcelain", "--untracked-files=no"):
        commit += " (+ local changes)"
    results = run_all(
        tests, env, Path(OUTDIR, "logs", "current"), args.timeout, args.break_on_errors
    )
    previous = load_json(CURRENT_JSON) if args.only else None
    save(CURRENT_JSON, branch, commit, tests, results, previous)
    if REFERENCE_JSON.exists():
        write_report(args)
    else:
        printf("No reference yet: run `init` to build it from main.")


if __name__ == "__main__":
    main()
