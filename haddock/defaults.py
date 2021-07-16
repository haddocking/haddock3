"""All default parameters used by the framework"""
import os
import multiprocessing
from pathlib import Path


# Locate the CNS binary
CNS_EXE = os.getenv("HADDOCK3_CNS_EXE")
if not CNS_EXE:
    bin_path = Path(__file__).resolve().parent.parent.absolute()
    CNS_EXE = bin_path / "bin/cns/cns_solve-1.31-UU-MacIntel.exe"

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

    data_path = Path(__file__).resolve().parent.absolute() / "data"

    PARAMETERS_FILE = data_path / "toppar/haddock.param"

    TOPOLOGY_FILE = data_path / "toppar/haddock.top"

    LINK_FILE = data_path / "toppar/protein-allhdg5-4-noter.link"

    TRANSLATION_VECTORS = {}
    for i in range(51):
        _s = f'trans_vector_{i}'
        _p = Path(data_path, 'toppar', 'initial_positions', _s)
        TRANSLATION_VECTORS[_s] = _p

    TENSORS = {
        "tensor_psf": data_path / "toppar/tensor.psf",
        "tensor_pdb": data_path / "toppar/tensor.pdb",
        "tensor_para_psf": data_path / "toppar/tensor_para.psf",
        "tensor_para_pdb": data_path / "toppar/tensor_para.pdb",
        "tensor_dani_psf": data_path / "toppar/tensor_dani.psf",
        "tensor_dani_pdb": data_path / "toppar/tensor_dani.pdb"
    }

    SCATTER_LIB = data_path / "toppar/scatter.lib"

    AXIS = {
        "top_axis": data_path / "toppar/top_axis.pro",
        "par_axis": data_path / "toppar/par_axis.pro",
        "top_axis_dani": data_path / "toppar/top_axis_dani.pro"
    }

    WATER_BOX = {
        "boxtyp20": data_path / "toppar/boxtyp20.pdb"
    }
