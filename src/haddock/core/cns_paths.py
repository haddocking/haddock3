"""Path to CNS-related files."""
from pathlib import Path

from haddock import toppar_path


PARAMETERS_FILE = "haddock.param"
TOPOLOGY_FILE = "haddock.top"
LINK_FILE = "protein-allhdg5-4-noter.link"
SCATTER_LIB = "scatter.lib"


def get_translation_vectors(path):
    """Generate paths for translation vectors."""
    translation_vectors = {}
    for i in range(51):
        _s = f'trans_vector_{i}'
        _p = Path(path, 'initial_positions', _s)
        translation_vectors[_s] = _p

    return translation_vectors


def get_tensors(path):
    """Generate paths for tensors."""
    tensors = {
        "tensor_psf": Path(path, "tensor.psf"),
        "tensor_pdb": Path(path, "tensor.pdb"),
        "tensor_para_psf": Path(path, "tensor_para.psf"),
        "tensor_para_pdb": Path(path, "tensor_para.pdb"),
        "tensor_dani_psf": Path(path, "tensor_dani.psf"),
        "tensor_dani_pdb": Path(path, "tensor_dani.pdb"),
        }
    return tensors


def get_axis(path):
    """Generate paths for axis."""
    axis = {
        "top_axis": Path(path, "top_axis.pro"),
        "par_axis": Path(path, "par_axis.pro"),
        "top_axis_dani": Path(path, "top_axis_dani.pro"),
        }
    return axis


def get_water_box(path):
    """Generate paths for water box."""
    water_box = {
        "boxtyp20": Path(path, "boxtyp20.pdb"),
        }
    return water_box


axis = get_axis(toppar_path)
link_file = Path(toppar_path, LINK_FILE)
parameters_file = Path(toppar_path, PARAMETERS_FILE)
scatter_lib = Path(toppar_path, SCATTER_LIB)
tensors = get_tensors(toppar_path)
topology_file = Path(toppar_path, TOPOLOGY_FILE)
translation_vectors = get_translation_vectors(toppar_path)
water_box = get_water_box(toppar_path)
