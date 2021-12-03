"""
Processes input PDB files to ensure compatibility with HADDOCK3.


"""
# search for ABSTRACT to reach the section where abstractions are defined
# search for CHECKING to reach the section where checking functions are defined
# search for CORRECT to reach the section where checking functions are defined
#
# functions can be defined with `def` or as variables with functools.partial
#
import itertools as it
import os
import string
from pathlib import Path
from functools import partial

from pdbtools import (
    pdb_selaltloc,
    pdb_rplresname,
    pdb_fixinsert,
    pdb_keepcoord,
    pdb_occ,
    pdb_tidy,
    pdb_segxchain,
    pdb_chainxseg,
    pdb_chain,
    pdb_reres,
    pdb_reatom,
    )

from haddock import log
from haddock.core.supported_molecules import (
    supported_cofactors,
    supported_carbohydrates,
    supported_ions,
    supported_modified_amino_acids,
    supported_multiatom_ions,
    supported_natural_amino_acids,
    supported_nucleic_acid_bases,
    )
from haddock.libs.libio import open_files_to_lines
from haddock.libs.libfunc import chainf
from haddock.libs.libpdb import read_chainids, read_segids, slc_resname


_OSUFFIX = '_processed'
_upper_case = list(string.ascii_uppercase)[::-1]
_CHAINS = it.cycle(_upper_case)


supported_hetatm = set(it.chain(
    supported_carbohydrates,
    supported_multiatom_ions,
    supported_cofactors,
    supported_ions,
    ))
supported_atom = set(it.chain(
    supported_modified_amino_acids,
    supported_natural_amino_acids,
    supported_nucleic_acid_bases,
    ))


def _open_or_give(lines_or_paths):
    """
    Adapt input to the functions.

    Parameters
    ----------
    lines_or_paths : list
        If is a list of strings (expected PDB lines), return it as is.
        If is a list of pathlib.Paths, opens them and returns a list
        of strings (lines).

    Raises
    ------
    TypeError
        In any other circumstances.
    """
    if all(isinstance(i, Path) for i in lines_or_paths):
        return open_files_to_lines(*lines_or_paths)
    elif all(isinstance(i, (list, tuple)) for i in lines_or_paths):
        return lines_or_paths
    else:
        raise TypeError('Unexpected types in `structures`')


def check_and_process_pdbs(paths_or_lines, **params):
    """Check and process PDBs according to HADDOCK3's specifications."""
    structures = _open_or_give(paths_or_lines)
    check_pdbs(structures, **params)
    result = process_pdbs(
        structures,
        save_output=params.pop('save_output', False),
        osuffix=params.pop('osuffix', _OSUFFIX),
        **params)
    return result


def check_pdbs(paths_or_lines, **param):
    """
    Check PDBs according to HADDOCK3 specifications.
    """
    structures = _open_or_give(paths_or_lines)

    # checks and log messages (do not alter the lines)
    # input should be list and not a generator
    individual_checks = [
        partial(check_supported_hetatm, user_defined=param),
        partial(check_supported_atom, user_defined=param),
        ]

    for structure, func in it.product(structures, individual_checks):
        func(structure)

    # these are the checks that are performed considering all the PDBs
    # combined
    combined_checks = [
        ]

    # perform the combined checking steps
    for check_func in combined_checks:
        check_func(structures)

    return


def process_pdbs(
        paths_or_lines,
        save_output=False,
        osuffix=_OSUFFIX,
        **param):
    """
    Processes PDB file contents for HADDOCK3 compatibility.

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
    structures = _open_or_give(paths_or_lines)

    # these are the processing or checking functions that should (if needed)
    # modify the input PDB and return the corrected lines
    # (in the like of pdbtools)
    line_by_line_processing_steps = [
        pdb_keepcoord.run,
        pdb_selaltloc.run,
        partial(pdb_occ.run, occupancy=1.00),
        replace_MSE_to_MET,
        replace_HSD_to_HIS,
        replace_HSE_to_HIS,
        replace_HID_to_HIS,
        replace_HIE_to_HIS,
        #partial(pdb_fixinsert.run, option_list=[]),
        ##
        partial(remove_unsupported_hetatm, user_defined=param),
        partial(remove_unsupported_atom, user_defined=param),
        ##
        partial(pdb_reatom.run, starting_value=1),
        partial(pdb_reres.run, starting_resid=1),
        #partial(pdb_tidy.run, strict=True),
        ##
        partial(map, lambda x: x.rstrip(os.linesep)),
        ]

    # perform the individual processing steps
    processed_individually = [
        list(chainf(structure, *line_by_line_processing_steps))
        for structure in structures
        ]

    whole_pdb_processing_steps = [
        #solve_no_chainID_no_segID,
        ]

    processed_combined_steps = [
        correct_equal_chain_segids,
        ]

    processed_combined = chainf(structures, *processed_combined_steps)

    if save_output:
        try:  # in case paths_or_lines is paths
            o_paths = (
                Path(f.parent, f.stem + osuffix).with_suffix(f.suffix)
                for f in map(Path, paths_or_lines)
                )
        except TypeError:
            o_paths = map(
                'processed_{}.pdb'.format,
                range(1, len(structures) + 1),
                )

        for out_p, lines in zip(o_paths, processed_combined):
            out_p.write_text(os.linesep.join(lines) + os.linesep)

        return

    else:
        return processed_individually


# CHECKING block


def check_supported_molecules(
        lines,
        haddock3_defined=None,
        user_defined=None,
        line_startswith=('ATOM', 'HETATM'),
        ):
    """."""
    user_defined = user_defined or set()
    haddock3_defined = haddock3_defined or set()
    allowed = set.union(haddock3_defined, user_defined)

    not_allowed_found = set()

    for line in lines:
        if line.startswith(line_startswith):
            residue = line[slc_resname].strip()
            if residue not in allowed:
                not_allowed_found.add(residue)

    for i in not_allowed_found:
        log.warning(
            f'* WARNING: residue {i!r} found and not supported as '
            f'{line_startswith!r}. This molecule will be removed.'
            )

    return


check_supported_hetatm = partial(
    check_supported_molecules,
    haddock3_defined=supported_hetatm,
    line_startswith='HETATM',
    )

check_supported_atom = partial(
    check_supported_molecules,
    haddock3_defined=supported_atom,
    line_startswith='ATOM',
    )


# CORRECT


def remove_unsupported_molecules(
        lines,
        haddock3_defined=None,
        user_defined=None,
        line_startswith=('ATOM', 'HETATM'),
        ):
    """."""
    user_defined = user_defined or set()
    haddock3_defined = haddock3_defined or set()
    allowed = set.union(haddock3_defined, user_defined)

    not_allowed_found = set()

    for line in lines:
        if line.startswith(line_startswith):
            residue = line[slc_resname].strip()
            if residue in allowed:
                yield line
            else:
                continue
        else:
            yield line

    return


remove_unsupported_hetatm = partial(
    remove_unsupported_molecules,
    haddock3_defined=supported_hetatm,
    line_startswith='HETATM',
    )

remove_unsupported_atom = partial(
    remove_unsupported_molecules,
    haddock3_defined=supported_atom,
    line_startswith='ATOM',
    )


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
        _chains = pdb_chain.run(lines, next(_CHAINS))
        new_lines = pdb_chainxseg.run(_chains)

    # gives priority to chains
    elif chainids != segids:
        new_lines = pdb_chainxseg.run(lines)

    else:
        # nothing to do
        return lines

    return list(new_lines)


def replace_HETATM_to_ATOM(fhandler, res):
    """."""
    for line in fhandler:
        if line.startswith('HETATM') and line[slc_resname].strip() == res:
            yield 'ATOM  ' + line[6:]
        else:
            yield line


def replace_residue(fhandler, resin, resout):
    """Replace residue by another and changes HETATM to ATOM if needed."""
    _ = replace_HETATM_to_ATOM(fhandler, res=resin)
    yield from pdb_rplresname.run(_, name_from=resin, name_to=resout)


replace_MSE_to_MET = partial(replace_residue, resin='MSE', resout='MET')
replace_HSD_to_HIS = partial(replace_residue, resin='HSD', resout='HIS')
replace_HSE_to_HIS = partial(replace_residue, resin='HSE', resout='HIS')
replace_HID_to_HIS = partial(replace_residue, resin='HID', resout='HIS')
replace_HIE_to_HIS = partial(replace_residue, resin='HIE', resout='HIS')


def correct_equal_chain_segids(structures):
    """
    Parameters
    ----------
    structures : list
    """


    return processed_structures
