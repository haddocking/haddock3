
"""
Library of functions to perform sequence and structural alignments.

Main functions
--------------

* :py:func:`write_unclustered_list`
"""
import os
from pathlib import Path
import time
from haddock import log


def write_unclustered_list(input_models, clustered_models):
    """
    Gets the list of unclustered structures.
    
    Parameters
    ----------
    input_models : list
        list of input models
    clustered_models : list
        list of clustered models
    """
    start = time.time()
    output_fname = Path('unclustered.txt')
    output_str = f'### Unclustered structures ###{os.linesep}'
    output_str += os.linesep
    unclustered_list = []
    # checking which input models have not been clustered
    for model in input_models:
        if model not in clustered_models:
            unclustered_list.append((model.score, model))
    unclustered_list.sort()
    # adding models to output string
    for model_ranking, sorted_uncl in enumerate(unclustered_list, start=1):
        output_str += f'{model_ranking}\t{sorted_uncl[1].file_name}\t{sorted_uncl[0]:.2f}{os.linesep}'
    output_str += os.linesep
    log.info('Saving unclustered structures to unclustered.txt')
    with open(output_fname, 'w') as out_fh:
            out_fh.write(output_str)
    elap_time = time.time() - start
    log.info(f"unclustered structures written in {elap_time} seconds")
    
