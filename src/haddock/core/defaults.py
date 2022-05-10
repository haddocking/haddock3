"""All default parameters used by the framework."""
import string
import sys
from pathlib import Path

import yaml

from haddock import core_path, haddock3_repository_path, log


# Locate the CNS binary
cns_exec = Path(haddock3_repository_path, "bin", "cns")
if not cns_exec.exists():
    log.error(
        'CNS executable `bin/cns` not found. '
        'Did you installed HADDOCK3 properly?'
        )
    sys.exit()


# Module input and generated data will be stored in folder starting by
#  this prefix
MODULE_PATH_NAME = "step_"

# Default name for exchange module information file
MODULE_IO_FILE = "io.json"

# Temptative number of max allowed number of modules to execute
MAX_NUM_MODULES = 10000

valid_run_dir_chars = string.ascii_letters + string.digits + "._-/\\"

RUNDIR = "run_dir"

with open(Path(core_path, "mandatory.yaml"), 'r') as fin:
    _ycfg = yaml.safe_load(fin)
max_molecules_allowed = _ycfg["molecules"]["maxitems"]
del _ycfg
