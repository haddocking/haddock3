"""Path to CNS-related files."""
from pathlib import Path

from haddock import toppar_path


parameters_file = Path(toppar_path, "haddock.param")

topology_file = Path(toppar_path, "haddock.top")

link_file = Path(toppar_path, "protein-allhdg5-4-noter.link")

translation_vectors = {}
for i in range(51):
    _s = f'trans_vector_{i}'
    _p = Path(toppar_path, 'initial_positions', _s)
    translation_vectors[_s] = _p

tensors = {
    "tensor_psf": Path(toppar_path, "tensor.psf"),
    "tensor_pdb": Path(toppar_path, "tensor.pdb"),
    "tensor_para_psf": Path(toppar_path, "tensor_para.psf"),
    "tensor_para_pdb": Path(toppar_path, "tensor_para.pdb"),
    "tensor_dani_psf": Path(toppar_path, "tensor_dani.psf"),
    "tensor_dani_pdb": Path(toppar_path, "tensor_dani.pdb"),
    }

scatter_lib = Path(toppar_path, "scatter.lib")

axis = {
    "top_axis": Path(toppar_path, "top_axis.pro"),
    "par_axis": Path(toppar_path, "par_axis.pro"),
    "top_axis_dani": Path(toppar_path, "top_axis_dani.pro"),
    }

water_box = {
    "boxtyp20": Path(toppar_path, "boxtyp20.pdb"),
    }
