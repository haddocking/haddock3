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
    pdb_segxchain,
    pdb_chainxseg,
    pdb_chain,
    )

from haddock import log
from haddock.core.supported_molecules import (
    supported_carbohydrates,
    supported_ions,
    supported_modified_amino_acids,
    supported_multiatom_ions,
    supported_nucleic_acids_bases,
    )
from haddock.libs.libio import open_files_to_lines
from haddock.libs.libfunc import chainf
from haddock.libs.libpdb import read_chainids, read_segids, slc_resname


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

    procfiles = process_pdbs(files_content)

    # prepares output paths
    o_paths = [
        Path(f.parent, f.stem + osuffix).with_suffix(f.suffix)
        for f in map(Path, files)
        ]

    # saves processed structures to files
    for out_p, lines in zip(o_paths, procfiles):
        out_p.write_text(''.join(lines))

    return



def process_pdbs(structures):
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
    # checks and log messages (do not alter the lines)
    # input should be list and not a generator
    checks = [
        partial(check_supported_carbohydrates, user_defined=param),
        partial(check_supported_ions, user_defined=param),
        partial(check_supported_nucleic_acids, user_defined=param),
        partial(check_supported_modified_aa, user_defined=param),
        ]

    for structure in structures:
        for func in checks:
            func(structure)



    # these are the processing or checking functions that should (if needed)
    # modify the input PDB and return the corrected lines
    # (in the like of pdbtools)
    line_by_line_processing_steps = [
        pdb_keepcoord.run,
        pdb_selaltloc.run,
        partial(pdb_rplresname.run, name_from='MSE', name_to='MET'),
        partial(pdb_fixinsert.run, option_list=[]),
        #
        partial(remove_unsupported_carbohydrates, user_defined=param),
        partial(remove_unsupported_ions, user_defined=param),
        partial(remove_unsupported_nucleic_acids, user_defined=param),
        partial(remove_unsupported_modified_aa, user_defined=param),
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


def check_supported_molecules(
        lines,
        haddock3_defined=None,
        user_defined=None,
        ):
    """."""
    user_defined = user_defined or tuple()
    haddock3_defined = haddock3_defined or tuple()
    allowed = set(haddock3_defined) + set(user_defined)

    not_allowed_found = set()

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            residue = line[slc_resname]
            if residue not in allowed:
                not_allowed_found.add(residue)

    for i in not_allowed_found:
        log.warning(
            f'* WARNING: {i!r} found and not supported. This '
            'molecule will be removed.'
            )

    return


check_supported_carbohydrates = partial(check_supported_molecules, haddock3_defined=supported_carbohydrates)  # noqa: E221,E501
check_supported_ions          = partial(check_supported_molecules, haddock3_defined=supported_ions + supported_multiatom_ions)  # noqa: E221,E501
check_supported_nucleic_acids = partial(check_supported_molecules, haddock3_defined=supported_nucleic_acids_bases)  # noqa: E221,E501
check_supported_modified_aa   = partial(check_supported_molecules, haddock3_defined=supported_modified_amino_acids)  # noqa: E221,E501


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
    elif chainids != segids:
        new_lines = pdb_chainxseg.run(lines)

    else:
        # nothing to do
        return lines

    return list(new_lines)
