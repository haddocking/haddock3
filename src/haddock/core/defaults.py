"""All default parameters used by the framework"""
import os
import sys
import multiprocessing
import logging
from pathlib import Path

from haddock import haddock3_source_path, toppar_path


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

# Topology got a first-class treatment concerning folder structure
TOPOLOGY_PATH = "topology"

# Temptative number of max allowed number of modules to execute
MAX_NUM_MODULES = 10000


class Default:

    PARAMETERS_FILE = Path(toppar_path, "haddock.param")

    TOPOLOGY_FILE = Path(toppar_path, "haddock.top")

    LINK_FILE = Path(toppar_path, "protein-allhdg5-4-noter.link")

    TRANSLATION_VECTORS = {}
    for i in range(51):
        _p = Path(topopar, 'initial_positions', f'trans_vector_{i}')
        TRANSLATION_VECTORS[_s] = _p

    TENSORS = {
        "tensor_psf": Path(toppar, "tensor.psf"),
        "tensor_pdb": Path(toppar, "tensor.pdb"),
        "tensor_para_psf": Path(toppar, "tensor_para.psf"),
        "tensor_para_pdb": Path(toppar, "tensor_para.pdb"),
        "tensor_dani_psf": Path(toppar, "tensor_dani.psf"),
        "tensor_dani_pdb": Path(toppar, "tensor_dani.pdb"),
        }

    SCATTER_LIB = Path(toppar, "scatter.lib")

    AXIS = {
        "top_axis": Path(toppar, "top_axis.pro"),
        "par_axis": Path(toppar, "par_axis.pro"),
        "top_axis_dani": Path(toppar, "top_axis_dani.pro"),
        }

    WATER_BOX = {
        "boxtyp20": Path(toppar, "boxtyp20.pdb")
        }
