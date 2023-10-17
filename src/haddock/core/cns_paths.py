"""
Path to CNS-related files.

Most paths are defined by dictionaries that gather several related
paths. Here, instead of defining the dictionaries with static paths, we
have functions that create those dict-containing paths dynamically. The
default values are defined by:

- axis
- tensors
- translation_vectors
- water_box

But you can re-use the functions to create new dictionaries with updated
paths. This is useful for those cases when the `cns/` folder is moved
to a different folder.
"""
from pathlib import Path

from haddock import toppar_path
from haddock.core.typing import FilePath


# exact file names as present in the cns/ scripts folder
LINK_FILE = "protein-allhdg5-4-noter.link"
SCATTER_LIB = "scatter.lib"
INITIAL_POSITIONS_DIR = "initial_positions"

# default prepared paths
link_file = Path(toppar_path, LINK_FILE)
scatter_lib = Path(toppar_path, SCATTER_LIB)


def get_translation_vectors(path: FilePath) -> dict[str, Path]:
    """
    Generate paths for translation vectors.

    Parameters
    ----------
    path : pathlib.Path
        If absolute, paths will be absolute, if relative paths will be
        relative. Adds the INITIAL_POSITIONS_DIR path before the file
        name.
    """
    translation_vectors: dict[str, Path] = {}
    for i in range(51):
        _s = f'trans_vector_{i}'
        _p = Path(path, INITIAL_POSITIONS_DIR, _s)
        translation_vectors[_s] = _p

    return translation_vectors


def get_tensors(path: FilePath) -> dict[str, Path]:
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


def get_axis(path: FilePath) -> dict[str, Path]:
    """Generate paths for axis."""
    axis = {
        "top_axis": Path(path, "top_axis.pro"),
        "par_axis": Path(path, "par_axis.pro"),
        "top_axis_dani": Path(path, "top_axis_dani.pro"),
        }
    return axis


def get_water_box(path: FilePath) -> dict[str, Path]:
    """Generate paths for water box."""
    water_box = {
        "boxtyp20": Path(path, "boxtyp20.pdb"),
        }
    return water_box


axis = get_axis(toppar_path)
tensors = get_tensors(toppar_path)
translation_vectors = get_translation_vectors(toppar_path)
water_box = get_water_box(toppar_path)
