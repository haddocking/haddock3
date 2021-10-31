"""All default parameters used by the framework"""
import os
import sys
import multiprocessing
import logging
from pathlib import Path

from haddock import haddock3_repository_path, haddock3_source_path


logger = logging.getLogger(__name__)

# Locate the CNS binary
cns_exec = Path(haddock3_repository_path, "bin", "cns")
if not cns_exec.exists():
    logger.error(
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
