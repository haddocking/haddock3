"""All default parameters used by the framework."""
import sys
from pathlib import Path

from haddock import haddock3_repository_path, log


# Locate the CNS binary
cns_exec = Path(haddock3_repository_path, "bin", "cns")
if not cns_exec.exists():
    log.error(
        'CNS executable `bin/cns` not found. '
        'Did you installed HADDOCK3 properly?'
        )
    sys.exit()

print(cns_exec)

# Module input and generated data will be stored in folder starting by
#  this prefix
MODULE_PATH_NAME = "step_"

# Default name for exchange module information file
MODULE_IO_FILE = "io.json"

# Temptative number of max allowed number of modules to execute
MAX_NUM_MODULES = 10000
