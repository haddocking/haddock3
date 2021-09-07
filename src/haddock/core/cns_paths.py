"""Path to CNS-related files."""
from pathlib import Path


parameters_file = Path(toppar_path, "haddock.param")

topology_file = Path(toppar_path, "haddock.top")

link_file = Path(toppar_path, "protein-allhdg5-4-noter.link")

translation_vectors = {}
for i in range(51):
    _s = f'trans_vector_{i}'
    _p = Path(topopar, 'initial_positions', _s)
    TRANSLATION_VECTORS[_s] = _p

tensors = {
    "tensor_psf": Path(toppar, "tensor.psf"),
    "tensor_pdb": Path(toppar, "tensor.pdb"),
    "tensor_para_psf": Path(toppar, "tensor_para.psf"),
    "tensor_para_pdb": Path(toppar, "tensor_para.pdb"),
    "tensor_dani_psf": Path(toppar, "tensor_dani.psf"),
    "tensor_dani_pdb": Path(toppar, "tensor_dani.pdb"),
    }

scatter_lib = Path(toppar, "scatter.lib")

axis = {
    "top_axis": Path(toppar, "top_axis.pro"),
    "par_axis": Path(toppar, "par_axis.pro"),
    "top_axis_dani": Path(toppar, "top_axis_dani.pro"),
    }

water_box = {
    "boxtyp20": Path(toppar, "boxtyp20.pdb")
    }
