"""Logic pertraining to the reading of the IO."""
from haddock.libs.libontology import Format
from haddock.libs.libontology import PDBFile


def load_from_previous(io_output, check_balance=False):
    """Wrapper function to get data from previous IO."""
    # Get the models generated in previous step
    models_to_refine = []
    if not all([isinstance(e, PDBFile) for e in io_output]):
        len_list = []
        for i, mol_dic in enumerate(io_output):
            len_list.append(len(mol_dic.values()))
            # sub_lists.append(list(mol_dic.values()))

        # check if it has the same dize
        if not all(x == len_list[0] for x in len_list):
            return False
        else:
            models_to_refine = []
            for element_key in io_output[0]:
                p = io_output[0][element_key]
                models_to_refine.append(p)
    else:
        models_to_refine = [
            p
            for p in io_output
            if p.file_type == Format.PDB
            ]

    return models_to_refine
