"""Gear to handle the reading of the IO."""
from haddock.libs.libontology import Format, PDBFile


def load_from_previous(io_output):
    """Load previous IO."""
    models_to_refine = []
    if not all(isinstance(e, PDBFile) for e in io_output):
        models_to_refine = []
        for mol_dic in io_output:
            for k in mol_dic:
                p = mol_dic[k]
                models_to_refine.append(p)

    else:
        models_to_refine = [
            p
            for p in io_output
            if p.file_type == Format.PDB
            ]

    return models_to_refine
