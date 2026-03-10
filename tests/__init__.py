"""Define common test variables."""

import shutil
import pytest
import platform
from pathlib import Path

from haddock.libs.libgrid import ping_dirac
from haddock.modules import modules_category
import platform as _platform

tests_path = Path(__file__).resolve().parents[0]
data_folder = Path(tests_path, "data")
golden_data = Path(tests_path, "golden_data")

# preprocessing files
broken_pdb = Path(data_folder, "broken.pdb")
corrected_pdb = Path(data_folder, "corrected.pdb")

residues_top = Path(data_folder, "residues.top")

configs_data = Path(tests_path, "configs")
emptycfg = Path(configs_data, "empty.cfg")
haddock3_yaml_cfg_examples = Path(configs_data, "yml_example.yml")
haddock3_yaml_converted = Path(configs_data, "yaml2cfg_converted.cfg")
haddock3_yaml_converted_no_header = Path(
    configs_data, "yaml2cfg_converted_no_header.cfg"
)
clean_steps_folder = Path(tests_path, "clean_output_data")

steptmp = Path(data_folder, "0_dummystep")

# defines which modules are already working
working_modules = [t for t in modules_category.items() if t[0] != "topocg"]


try:
    import py3Dmol

    NOTEBOOK_ENABLED = True
except ImportError:
    NOTEBOOK_ENABLED = False

has_notebook = pytest.mark.skipif(
    not NOTEBOOK_ENABLED, reason="notebook dependencies not found"
)

has_grid = pytest.mark.skipif(not ping_dirac(), reason="Dirac not reachable")

try:
    import deeprank_gnn.predict  # noqa: F401

    DEEPRANK_ENABLED = True
except (ImportError, ModuleNotFoundError):
    DEEPRANK_ENABLED = False

has_deeprank = pytest.mark.skipif(
    not DEEPRANK_ENABLED,
    reason="deeprank_gnn is not installed or not supported on this platform",
)

_CHROME_BINS = (
    "google-chrome",
    "google-chrome-stable",
    "chromium-browser",
    "chromium",
    "chrome",
)
has_chrome = pytest.mark.skipif(
    not any(shutil.which(b) for b in _CHROME_BINS),
    reason="Google Chrome not found (required by Kaleido for PNG export)",
)
is_linux_x86_64 = pytest.mark.skipif(
    platform.system().lower() != "linux" or platform.machine().lower() != "x86_64",
    reason="Only runs on x86_64 Linux systems",
)
