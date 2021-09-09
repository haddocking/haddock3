"""All default parameters used by the framework"""
import os
import sys
import multiprocessing
import logging
from pathlib import Path


logger = logging.getLogger(__name__)

# Locate the CNS binary
CNS_EXE = os.getenv("HADDOCK3_CNS_EXE")
if not CNS_EXE:
    bin_path = Path(__file__).resolve().parent.parent.parent.absolute()
    CNS_EXE = bin_path / "bin" / "cns"
    if not CNS_EXE.exists():
        logger.error('HADDOCK3_CNS_EXE not defined and bin/cns not found')
        sys.exit()

# Number of cores to use
NUM_CORES = int(os.getenv("HADDOCK3_NUM_CORES", multiprocessing.cpu_count()))

# Module input and generated data will be stored in folder starting by
#  this prefix
MODULE_PATH_NAME = "step_"

# Default name for exchange module information file
MODULE_IO_FILE = "io.json"

# Temptative number of max allowed number of modules to execute
MAX_NUM_MODULES = 10000
