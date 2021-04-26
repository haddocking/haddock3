import os
import multiprocessing
from pathlib import Path


# Locate the CNS binary
CNS_EXE = os.getenv('HADDOCK3_CNS_EXE')
if not CNS_EXE:
    bin_path = Path(__file__).resolve().parent.parent.absolute()
    CNS_EXE = bin_path / "bin/cns/cns_solve-1.31-UU-MacIntel.exe"

# Number of cores to use
try:
    NUM_CORES = int(os.getenv('HADDOCK3_NUM_CORES'))
    if not NUM_CORES:
        NUM_CORES = multiprocessing.cpu_count()
except ValueError:
    NUM_CORES = multiprocessing.cpu_count()

# Module input and generated data will be stored in folder starting by this prefix
MODULE_PATH_NAME = "step_"

# Default name for exchange module information file
MODULE_IO_FILE = "io.json"

# Temptative number of max allowed number of modules to execute
MAX_NUM_MODULES = 10000


class Default:

    data_path = Path(__file__).resolve().parent.absolute() / 'data'

    PARAMETERS_FILE = data_path / "toppar/haddock.param"

    TOPOLOGY_FILE = data_path / "toppar/haddock.top"

    LINK_FILE = data_path / "toppar/protein-allhdg5-4-noter.link"

    TRANSLATION_VECTORS = {
        "trans_vector_0": data_path / "toppar/initial_positions/trans_vector_0",
        "trans_vector_1": data_path / "toppar/initial_positions/trans_vector_1",
        "trans_vector_2": data_path / "toppar/initial_positions/trans_vector_2",
        "trans_vector_3": data_path / "toppar/initial_positions/trans_vector_3",
        "trans_vector_4": data_path / "toppar/initial_positions/trans_vector_4",
        "trans_vector_5": data_path / "toppar/initial_positions/trans_vector_5",
        "trans_vector_6": data_path / "toppar/initial_positions/trans_vector_6",
        "trans_vector_7": data_path / "toppar/initial_positions/trans_vector_7",
        "trans_vector_8": data_path / "toppar/initial_positions/trans_vector_8",
        "trans_vector_9": data_path / "toppar/initial_positions/trans_vector_9",
        "trans_vector_10": data_path / "toppar/initial_positions/trans_vector_10",
        "trans_vector_11": data_path / "toppar/initial_positions/trans_vector_11",
        "trans_vector_12": data_path / "toppar/initial_positions/trans_vector_12",
        "trans_vector_13": data_path / "toppar/initial_positions/trans_vector_13",
        "trans_vector_14": data_path / "toppar/initial_positions/trans_vector_14",
        "trans_vector_15": data_path / "toppar/initial_positions/trans_vector_15",
        "trans_vector_16": data_path / "toppar/initial_positions/trans_vector_16",
        "trans_vector_17": data_path / "toppar/initial_positions/trans_vector_17",
        "trans_vector_18": data_path / "toppar/initial_positions/trans_vector_18",
        "trans_vector_19": data_path / "toppar/initial_positions/trans_vector_19",
        "trans_vector_20": data_path / "toppar/initial_positions/trans_vector_20",
        "trans_vector_21": data_path / "toppar/initial_positions/trans_vector_21",
        "trans_vector_22": data_path / "toppar/initial_positions/trans_vector_22",
        "trans_vector_23": data_path / "toppar/initial_positions/trans_vector_23",
        "trans_vector_24": data_path / "toppar/initial_positions/trans_vector_24",
        "trans_vector_25": data_path / "toppar/initial_positions/trans_vector_25",
        "trans_vector_26": data_path / "toppar/initial_positions/trans_vector_26",
        "trans_vector_27": data_path / "toppar/initial_positions/trans_vector_27",
        "trans_vector_28": data_path / "toppar/initial_positions/trans_vector_28",
        "trans_vector_29": data_path / "toppar/initial_positions/trans_vector_29",
        "trans_vector_30": data_path / "toppar/initial_positions/trans_vector_30",
        "trans_vector_31": data_path / "toppar/initial_positions/trans_vector_31",
        "trans_vector_32": data_path / "toppar/initial_positions/trans_vector_32",
        "trans_vector_33": data_path / "toppar/initial_positions/trans_vector_33",
        "trans_vector_34": data_path / "toppar/initial_positions/trans_vector_34",
        "trans_vector_35": data_path / "toppar/initial_positions/trans_vector_35",
        "trans_vector_36": data_path / "toppar/initial_positions/trans_vector_36",
        "trans_vector_37": data_path / "toppar/initial_positions/trans_vector_37",
        "trans_vector_38": data_path / "toppar/initial_positions/trans_vector_38",
        "trans_vector_39": data_path / "toppar/initial_positions/trans_vector_39",
        "trans_vector_40": data_path / "toppar/initial_positions/trans_vector_40",
        "trans_vector_41": data_path / "toppar/initial_positions/trans_vector_41",
        "trans_vector_42": data_path / "toppar/initial_positions/trans_vector_42",
        "trans_vector_43": data_path / "toppar/initial_positions/trans_vector_43",
        "trans_vector_44": data_path / "toppar/initial_positions/trans_vector_44",
        "trans_vector_45": data_path / "toppar/initial_positions/trans_vector_45",
        "trans_vector_46": data_path / "toppar/initial_positions/trans_vector_46",
        "trans_vector_47": data_path / "toppar/initial_positions/trans_vector_47",
        "trans_vector_48": data_path / "toppar/initial_positions/trans_vector_48",
        "trans_vector_49": data_path / "toppar/initial_positions/trans_vector_49",
        "trans_vector_50": data_path / "toppar/initial_positions/trans_vector_50"
    }

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
