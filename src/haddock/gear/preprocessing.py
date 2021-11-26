"""
Processes input PDB files to ensure compatibility with HADDOCK3.


"""
import os
import string
from pathlib import Path
from functools import partial

from pdbtools import (
    pdb_selaltloc,
    pdb_rplresname,
    pdb_fixinsert,
    pdb_keepcoord,
    pdb_tidy,
    )

from haddock.libs.libfunc import chainf, chainfs, reduce_helper


upper_case = list(string.ascii_uppercase)[::-1]


def process_pdb_files(*files, osuffix='_processed'):
    """
    Processes and checks structures for HADDOCK3 compatibility.

    Parameters
    ----------
    files : str or pathlib.Path (many)
        Paths to structure files.

    osuffix : str
        A suffix to add to the new processed files when saved to disk.
        New files are saved in the same folder path as the input.

    Returns
    -------
    list of lists of strings
        The results from :func:`process_pdb`.
    """
    files_content = open_files_to_lines(*files)

    procfiles = process_structures(files_content)

    # prepares output paths
    f_paths = map(Path, files)
    o_paths = [
        Path(f.parent, f.stem + osuffix, f.suffix)
        for f in map(Path, files)
        ]

    # saves processed structures to files
    for out_p, lines in zip(o_paths, procfiles):
        out_p.write_text(''.join(lines))

    return



def process_pdb(structures):
    """
    Processes and checks PDB file contents for HADDOCK3 compatibility.

    Parameters
    ----------
    structures : list of lists of strings
        A list of lists, where each sublist corresponds to the lines of
        the structures (file content). Why paths are not accepted? So
        that `process_structure` can be used from loaded information.
        Why don't accept generators? Because we will need to sort at
        some point.

    Returns
    -------
    list of list
        The same data structure as `structures` but with the corrected
        (processed) content.
    """
    # these are the processing or checking functions that should (if needed)
    # modify the input PDB and return the corrected lines
    # (in the like of pdbtools)
    line_by_line_processing_steps = [
        pdb_keepcoord.run,
        pdb_selaltloc.run,
        partial(pdb_rplresname.run, name_from='MSE', name_to='MET'),
        partial(pdb_fixinsert.run, option_list=[]),
        #
        pdb_tidy.run,
        #
        #partial(map, vartial(str.rstrip, os.linesep)),
        partial(map, lambda x: x.rstrip(os.linesep)),
        ]

    whole_pdb_processing_steps = [
        solve_no_chainID_no_segID,
        ]

    # these are the checks that are performed considering all the PDBs
    # combined
    combined_checking_steps = [
        ]

    # perform the individual processing steps
    processed_individually = [
        list(chainf(structure, *line_by_line_processing_steps))
        for structure in structures
        ]

    # perform the combined checking steps
    for check_func in combined_checking_steps:
        check_func(processed_individually)

    return processed_individually


def solve_no_chainID_no_segID(lines):
    """
    Solve inconsistencies with chainID and segID.

    If segID is non-existant, copy chainID over segID, and vice-verse.
    If none are present, adds an upper case char starting from A. This
    char is not repeated until the alphabet exhausts.
    If chainIDs and segIDs differ, copy chainIDs over segIDs.

    Parameters
    ----------
    lines : list of str
        The lines of a PDB file.

    Returns
    -------
    list
        With new lines. Or the input ones if no modification was made.
    """
    chainids = read_chainids(lines)
    segids = read_segids(lines)

    if not chainids and segids:
        new_lines = pdb_segxchain.run(lines)

    elif chainids and not segids:
        new_lines = pdb_chainxseg.run(lines)

    elif not chainids and not segids:
        _chains = pdb_chain.run(lines, upper_case.pop())
        new_lines = pdb_chainxseg.run(_chains)

    # gives priority to chains
    elif chainds != segids:
        new_lines = pdb_chainxseg.run(lines)

    else:
        # nothing to do
        return lines

    return list(new_lines)
