"""
Library of functions related to the clustering modules.

Main functions
--------------

* :py:func:`write_unclustered_list`
"""
import os
from pathlib import Path

from haddock import log
from haddock.core.typing import FilePath
from haddock.libs.libontology import PDBFile


def write_structure_list(input_models: list[PDBFile],
                         clustered_models: list[PDBFile],
                         out_fname: FilePath) -> None:
    """
    Get the list of unclustered structures.
    
    Parameters
    ----------
    input_models : list
        list of input models
    clustered_models : list
        list of clustered models
    """
    output_fname = Path(out_fname)
    output_str = f'rank\tmodel_name\tscore\tcluster_id{os.linesep}'
    structure_list: list[PDBFile] = []
    # checking which input models have not been clustered
    for model in input_models:
        if model not in clustered_models:
            model.clt_id = "-"
            structure_list.append(model)
    # extending and sorting
    structure_list.extend(clustered_models)
    structure_list.sort(key=lambda model: model.score)
    # adding models to output string
    for mdl_rank, mdl in enumerate(structure_list, start=1):
        output_str += (
            f'{mdl_rank}\t{mdl.file_name}\t{mdl.score:.2f}\t{mdl.clt_id}'
            f'{os.linesep}'
            )
    output_str += os.linesep
    log.info(f'Saving structure list to {out_fname}')
    with open(output_fname, 'w') as out_fh:
        out_fh.write(output_str)
