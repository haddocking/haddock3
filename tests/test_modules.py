import importlib
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from haddock.libs.libgrid import GRIDScheduler
from haddock.libs.libhpc import HPCScheduler
from haddock.libs.libmpi import MPIScheduler
from haddock.libs.libparallel import Scheduler
from haddock.modules import get_engine


def test_get_engine_grid_unavailable():
    with patch("haddock.modules.ping_dirac", return_value=False):
        engine = get_engine(mode="grid", params={"ncores": 1, "max_cpus": 4})
        assert isinstance(engine, partial)
        assert engine.func is Scheduler


def test_get_engine_grid_available():
    with patch("haddock.modules.ping_dirac", return_value=True):
        engine = get_engine(mode="grid", params={"ncores": 1, "max_cpus": 4})
        assert isinstance(engine, partial)
        assert engine.func is GRIDScheduler


def test_get_engine_local():
    engine = get_engine(mode="local", params={"ncores": 2, "max_cpus": 8})
    assert isinstance(engine, partial)
    assert engine.func is Scheduler
    assert engine.keywords["ncores"] == 2
    assert engine.keywords["max_cpus"] == 8


def test_get_engine_batch():
    engine = get_engine(
        mode="batch",
        params={"queue": "short", "queue_limit": 100, "concat": 5},
    )
    assert isinstance(engine, partial)
    assert engine.func is HPCScheduler
    assert engine.keywords["target_queue"] == "short"
    assert engine.keywords["queue_limit"] == 100
    assert engine.keywords["concat"] == 5


def test_get_engine_mpi():
    engine = get_engine(mode="mpi", params={"ncores": 4})
    assert isinstance(engine, partial)
    assert engine.func is MPIScheduler
    assert engine.keywords["ncores"] == 4


@pytest.mark.parametrize(
    "module_name",
    [
        "haddock.modules.refinement.openmm",
        "haddock.modules.scoring.prodigyligand",
        "haddock.modules.scoring.prodigyprotein",
    ],
)
def test_non_cns_modules_run_locally_in_grid_mode(module_name, monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    module = importlib.import_module(module_name)
    step = module.HaddockModule(order=1, path=Path("."))
    step.params["mode"] = "grid"
    step.previous_io = SimpleNamespace(retrieve_models=lambda **_: [])
    with (
        patch.object(module, "get_engine", side_effect=RuntimeError) as get_engine,
        patch.object(step, "create_directories", create=True),
        patch.object(step, "set_max_cpu", create=True),
        pytest.raises(RuntimeError),
    ):
        step._run()
    assert get_engine.call_args.args[0] == "local"
