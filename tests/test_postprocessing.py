"""Test gear.postprocessing."""

import os
import tarfile

from haddock.gear import postprocessing
from haddock.gear.postprocessing import _robust_rmtree, archive_run


def _make_run_dir(base, name="run1"):
    run_dir = base / name
    run_dir.mkdir()
    (run_dir / "log").write_text("run log\n")
    (run_dir / "params.cfg").write_text("dummy\n")
    return run_dir


def test_archive_run_creates_archive_and_deletes(tmp_path, monkeypatch):
    """archive_run tarballs the run dir and removes the original."""
    monkeypatch.chdir(tmp_path)
    run_dir = _make_run_dir(tmp_path)

    run_archive, analysis_archive = archive_run(str(run_dir))

    assert os.path.exists(run_archive)
    assert tarfile.is_tarfile(run_archive)
    assert analysis_archive is None
    # original directory removed
    assert not run_dir.exists()


def test_archive_run_no_delete(tmp_path, monkeypatch):
    """archive_run keeps the run dir when delete=False."""
    monkeypatch.chdir(tmp_path)
    run_dir = _make_run_dir(tmp_path)

    archive_run(str(run_dir), delete=False)

    assert run_dir.exists()


def test_archive_run_survives_undeletable_dir(tmp_path, monkeypatch):
    """archive_run must not crash if the run dir cannot be fully removed."""
    monkeypatch.chdir(tmp_path)
    run_dir = _make_run_dir(tmp_path)

    # Simulate a non-local filesystem that never lets the tree be removed.
    def _always_enotempty(_):
        raise OSError(39, "Directory not empty")

    monkeypatch.setattr(postprocessing.shutil, "rmtree", _always_enotempty)

    run_archive, _ = archive_run(str(run_dir))
    assert os.path.exists(run_archive)


def test_robust_rmtree_removes_dir(tmp_path):
    """_robust_rmtree deletes a normal directory."""
    target = tmp_path / "gone"
    target.mkdir()
    (target / "f").write_text("x")

    _robust_rmtree(str(target))

    assert not target.exists()


def test_robust_rmtree_retries_then_gives_up(tmp_path, monkeypatch, caplog):
    """_robust_rmtree retries and logs a warning instead of raising."""
    calls = []

    def _always_enotempty(_):
        calls.append(1)
        raise OSError(39, "Directory not empty")

    monkeypatch.setattr(postprocessing.shutil, "rmtree", _always_enotempty)
    monkeypatch.setattr(postprocessing.time, "sleep", lambda _s: None)

    # Should not raise.
    _robust_rmtree("/some/run1", retries=3, delay=0)

    assert len(calls) == 3


def test_robust_rmtree_recovers_on_retry(tmp_path, monkeypatch):
    """_robust_rmtree succeeds if a later attempt works (transient failure)."""
    real_rmtree = postprocessing.shutil.rmtree
    target = tmp_path / "flaky"
    target.mkdir()
    state = {"n": 0}

    def _flaky(path):
        state["n"] += 1
        if state["n"] < 2:
            raise OSError(39, "Directory not empty")
        real_rmtree(path)

    monkeypatch.setattr(postprocessing.shutil, "rmtree", _flaky)
    monkeypatch.setattr(postprocessing.time, "sleep", lambda _s: None)

    _robust_rmtree(str(target))

    assert not target.exists()
    assert state["n"] == 2
