"""All default parameters used by the framework."""
import string
import sys
from pathlib import Path
import os

import yaml

from haddock import core_path, haddock3_repository_path, log


# Locate the CNS binary
cns_exec = os.getenv("CNS_EXEC")
if cns_exec is None:
    log.error(
        'CNS_EXEC not defined. '
        'Please check the install instructions.'
        )
    sys.exit()
cns_exec = Path(cns_exec)

if not cns_exec.exists():
    log.error(
        'CNS_EXEC points to a non-existing file. '
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
