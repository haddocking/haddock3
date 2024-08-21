from pathlib import Path

import pytest
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces


def dockq_wrapper(model: Path, native: Path) -> dict[str, float]:
    """Wrapper to run DockQ analysis on a model and a native structure."""
    native_structure = load_PDB(str(model))
    model_structure = load_PDB(str(native))
    dockq_analysis = run_on_all_native_interfaces(
        model_structure=model_structure, native_structure=native_structure
    )[0]

    return dockq_analysis["AB"]


@pytest.fixture
def run_dockq_analysis():
    return dockq_wrapper
