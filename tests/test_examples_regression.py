"""Test the example regression script `examples/run_examples_regression.py`."""

import copy
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(
    Path(__file__).resolve().parents[1],
    "examples",
    "run_examples_regression.py",
)

CAPRI_HEADER = ["model", "md5", "caprieval_rank", "score", "irmsd", "fnat",
                "lrmsd", "ilrmsd", "dockq"]

# rank, irmsd, fnat, ilrmsd, dockq
CAPRI_ROWS = [
    (1, 2.207, 0.389, 4.677, 0.473),  # acceptable
    (2, 1.372, 0.417, 2.653, 0.621),  # medium
    (3, 9.000, 0.000, 0.900, 0.050),  # incorrect (but ligand-medium)
    (4, 0.992, 0.667, 1.700, 0.812),  # high, lowest i-RMSD
    (5, "nan", "nan", "nan", "nan"),  # no reference values
]


@pytest.fixture(name="regression", scope="module")
def fixture_regression():
    """Load the regression script as a module."""
    spec = importlib.util.spec_from_file_location("regression", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_capri(path):
    """Write a small capri_ss.tsv file."""
    lines = ["\t".join(CAPRI_HEADER)]
    for rank, irmsd, fnat, ilrmsd, dockq in CAPRI_ROWS:
        values = [f"../model_{rank}.pdb", "-", rank, -100, irmsd, fnat, 0,
                  ilrmsd, dockq]
        lines.append("\t".join(map(str, values)))
    path.write_text("\n".join(lines) + "\n")


def make_run(run_dir, finished=True, error=False):
    """Create a minimal finished run directory."""
    for step in ("0_topoaa", "1_caprieval"):
        Path(run_dir, step).mkdir(parents=True)
        Path(run_dir, step, "io.json").write_text('{"input": [], "output": [{}, {}]}')  # noqa: E501
    write_capri(Path(run_dir, "1_caprieval", "capri_ss.tsv"))
    log = ["[2026-01-01 10:00:00,000 libworkflow INFO] Reading instructions"]
    if error:
        log.append("[2026-01-01 10:00:01,000 libutil ERROR] Something broke")
    if finished:
        log.append("[2026-01-01 10:00:02,000 cli INFO] Finished at 01/01/2026.")
    Path(run_dir, "log").write_text("\n".join(log) + "\n")


def test_discovery_skips(regression):
    """MPI, batch and follow-up examples are skipped, others are run."""
    tests = {t[0]: t[2] for t in regression.discover_tests(regression.EXAMPLES)}
    assert tests["docking-protein-protein/docking-protein-protein-test"] is None
    assert tests["docking-protein-protein/docking-protein-protein-hpc-test"] == "mode = batch"  # noqa: E501
    assert tests["docking-protein-protein/docking-restart-exit-test"]
    assert all(reason for tid, reason in tests.items() if "mpi" in tid)


def test_expected_steps_stop_at_exit(regression):
    """Steps after an [exit] module are not expected."""
    cfg = Path(regression.EXAMPLES, "docking-protein-protein", "docking-exit-test.cfg")  # noqa: E501
    assert regression.expected_steps(cfg) == ["topoaa", "rigidbody"]


def test_capri_stats_dockq(regression, tmp_path):
    """DockQ criteria statistics follow extract-capri-stats-dockq.sh."""
    capri = Path(tmp_path, "capri_ss.tsv")
    write_capri(capri)
    stats = regression.capri_stats(capri, regression.DOCKQ_CRITERIA)
    assert (stats["nmodels"], stats["acceptable"], stats["medium"], stats["high"]) == (5, 3, 2, 1)  # noqa: E501
    assert stats["first_acceptable"]["rank"] == "1"
    assert stats["first_medium"]["rank"] == "2"
    assert stats["best"]["rank"] == "4"
    assert stats["best"]["dockq"] == pytest.approx(0.812)


def test_capri_stats_ligand(regression, tmp_path):
    """Ligand criteria use il-RMSD and pick the best model by il-RMSD."""
    capri = Path(tmp_path, "capri_ss.tsv")
    write_capri(capri)
    stats = regression.capri_stats(capri, regression.LIGAND_CRITERIA)
    assert (stats["acceptable"], stats["medium"], stats["high"]) == (2, 1, 0)
    assert stats["first_acceptable"]["rank"] == "3"
    assert stats["best"]["rank"] == "3"


@pytest.mark.parametrize(
    "finished,error,expected",
    [
        (True, False, []),
        (False, False, ["run did not finish"]),
        (True, True, ["ERROR/CRITICAL"]),
    ],
)
def test_inspect_run(regression, tmp_path, finished, error, expected):
    """The completion check reports unfinished runs and log errors."""
    run_dir = Path(tmp_path, "run1")
    make_run(run_dir, finished=finished, error=error)
    result = regression.inspect_run(run_dir, ["topoaa", "caprieval"], "x")
    assert len(result["problems"]) == len(expected)
    for problem, text in zip(result["problems"], expected):
        assert text in problem
    assert result["steps"]["0_topoaa"]["nmodels"] == 2
    assert "dockq" in result["steps"]["1_caprieval"]["capri"]


def test_inspect_run_missing_step(regression, tmp_path):
    """A module without a step folder is reported."""
    run_dir = Path(tmp_path, "run1")
    make_run(run_dir)
    result = regression.inspect_run(run_dir, ["topoaa", "caprieval", "emref"], "x")  # noqa: E501
    assert result["problems"] == ["step 2_emref missing"]


def _results(regression, tmp_path):
    run_dir = Path(tmp_path, "run1")
    make_run(run_dir)
    test = regression.inspect_run(run_dir, ["topoaa", "caprieval"], "x")
    test.update(status="OK", runtime=100.0, log_tail=[], returncode=0)
    return {"tests": {"x/test": test}, "skipped": {}}


def test_compare_unchanged(regression, tmp_path):
    """Identical results are reported as unchanged."""
    ref = _results(regression, tmp_path)
    _, summary = regression.compare(ref, copy.deepcopy(ref), 0.05, 0.5, 0.5)
    assert summary == [("x/test", "OK")]


def test_compare_detects_changes(regression, tmp_path):
    """Quality, ranking and output changes are reported."""
    ref = _results(regression, tmp_path)
    cur = copy.deepcopy(ref)
    stats = cur["tests"]["x/test"]["steps"]["1_caprieval"]["capri"]["dockq"]
    stats["acceptable"] = 2
    stats["first_acceptable"]["rank"] = "3"
    stats["best"].update(rank="1", dockq=0.5)
    cur["tests"]["x/test"]["steps"]["0_topoaa"]["nmodels"] = 1

    report, summary = regression.compare(ref, cur, 0.05, 0.5, 0.5)
    text = "\n".join(report)
    assert summary == [("x/test", "CHANGED")]
    assert "acceptable: 3 → 2 ▼" in text
    assert "first acceptable rank: 1 → 3 ▼" in text
    assert "best model rank: 4 → 1 ▲" in text
    assert "1 output models (reference 2)" in text
    assert "rank of the first acceptable or better model: 0 better, 0 unchanged, 1 worse" in text  # noqa: E501


def test_compare_failures(regression, tmp_path):
    """New failures and failures on both sides are told apart."""
    ref = _results(regression, tmp_path)
    cur = copy.deepcopy(ref)
    cur["tests"]["x/test"].update(status="FAILED", problems=["exit code 1"])
    _, summary = regression.compare(ref, cur, 0.05, 0.5, 0.5)
    assert summary == [("x/test", "FAILED")]

    ref["tests"]["x/test"].update(status="FAILED", problems=["exit code 1"])
    _, summary = regression.compare(ref, cur, 0.05, 0.5, 0.5)
    assert summary == [("x/test", "FAILED (also on main)")]
