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
import string
from collections import namedtuple
from functools import partial, wraps
from pathlib import Path
from os import linesep

from pdbtools import (
    pdb_chain,
    pdb_chainxseg,
    pdb_fixinsert,
    pdb_keepcoord,
    pdb_occ,
    pdb_reatom,
    pdb_reres,
    pdb_rplresname,
    pdb_segxchain,
    pdb_selaltloc,
    pdb_tidy,
    )

from haddock import log
from haddock.core.supported_molecules import (
    read_supported_residues_from_top_file,
    supported_ATOM,
    supported_HETATM,
    supported_ions_atoms,
    supported_ions_elements,
    supported_ions_resnames,
    )
from haddock.libs.libio import open_files_to_lines
from haddock.libs.libfunc import chainf
from haddock.libs.libpdb import (
    read_chainids,
    read_segids,
    slc_resname,
    slc_model,
    slc_element,
    slc_name,
    )


_OSUFFIX = '_processed'
_upper_case = list(string.ascii_uppercase)[::-1]
_CHAINS = it.cycle(_upper_case)


def _allow_dry(log_msg):
    def decorator(function):
        @wraps(function)
        def wrapper(lines, *args, dry=False, **kwargs):

            if dry:
                in_lines = list(lines)
                result = list(function(in_lines, *args, **kwargs))

                # Here we could use sets to increase speed, but for the size
                # of the systems, we can actually use lists and get a sorted
                # result by default. I tried using difflib from STD. But it is
                # just too slow.

                # _ is line
                additions = [_ for _ in result if _ not in in_lines]
                deletions = [_ for _ in in_lines if _ not in result]

                la = len(additions)
                ld = len(deletions)

                add_lines = linesep.join(f'+ {_}' for _ in additions)
                del_lines = linesep.join(f'- {_}' for _ in deletions)

                log_msg_ = log_msg.format(*args, *kwargs.values())
                extended_log = (
                    f'[{log_msg_}] + {la} - {ld} lines',
                    add_lines,
                    del_lines,
                    )

                log.info(linesep.join(extended_log))
                return in_lines

            # on dry, maintain the generator functionality
            else:
                return function(lines, *args, **kwargs)

        return wrapper
    return decorator


# make pdb-tools dryrunable
wdry_pdb_selaltloc = _allow_dry('pdb_selaltloc')(pdb_selaltloc.run)
wdry_pdb_rplresname = _allow_dry('pdb_rplresname')(pdb_rplresname.run)
wdry_pdb_fixinsert = _allow_dry('pdb_fixinsert')(pdb_fixinsert.run)
wdry_pdb_keepcoord = _allow_dry('pdb_keepcoord')(pdb_keepcoord.run)
wdry_pdb_occ = _allow_dry('pdb_occ')(pdb_occ.run)
wdry_pdb_tidy = _allow_dry('pdb_tidy')(pdb_tidy.run)
wdry_pdb_segxchain = _allow_dry('pdb_segxchain')(pdb_segxchain.run)
wdry_pdb_chainxseg = _allow_dry('pbd_segxchain')(pdb_chainxseg.run)
wdry_pdb_chain = _allow_dry('pdb_chain')(pdb_chain.run)
wdry_pdb_reres = _allow_dry('pdb_reres')(pdb_reres.run)
wdry_pdb_reatom = _allow_dry('pdb_reatom')(pdb_reatom.run)
wdry_rstrip = _allow_dry("str.rstrip")(partial(map, lambda x: x.rstrip(linesep)))  # noqa: E501


def _open_or_give(paths_or_lines):
    """
    Adapt input to the functions.

    Parameters
    ----------
    paths_or_lines : list
        If is a list of strings (expected PDB lines), return it as is.
        If is a list of pathlib.Paths, opens them and returns a list
        of strings (lines).

    Raises
    ------
    TypeError
        In any other circumstances.
    """
    if all(isinstance(i, Path) for i in paths_or_lines):
        return open_files_to_lines(*paths_or_lines)
    elif all(isinstance(i, (list, tuple)) for i in paths_or_lines):
        return paths_or_lines
    else:
        raise TypeError('Unexpected types in `paths_or_lines`.')


def read_additional_residues(top_fname):
    """Read additional residues listed in a .top filename."""
    residues = read_supported_residues_from_top_file(
        top_fname,
        regex=r'\nRESIdue ([A-Z0-9]{1,3}).*\n',
        Residue=namedtuple('New', ['resname']),
        )

    return tuple(i.resname for i in residues)



def process_pdbs(
        paths_or_lines,
        save_output=False,
        osuffix=_OSUFFIX,
        dry=False,
        user_supported_residues=None,
        **param):
    """
    Process PDB file contents for HADDOCK3 compatibility.

    Parameters
    ----------
    paths_or_lines : list
        A list of paths pointing to files or a list of lists of strings
        with the lines of the file contents.

    save_output : bool
        Whether to directly save the output to files. If False, returns
        a list of lists with the altered contents.

    dry : bool
        Perform a dry run. That is, does not change anything, and just
        report.

    user_supported_residues : list, tuple, or set
        The new residues that are allowed.

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
        wdry_pdb_keepcoord,
        # wdry_pdb_selaltloc,
        partial(wdry_pdb_occ, occupancy=1.00),
        replace_MSE_to_MET,
        replace_HSD_to_HIS,
        replace_HSE_to_HIS,
        replace_HID_to_HIS,
        replace_HIE_to_HIS,
        add_charges_to_ions,
        convert_ATOM_to_HETATM,
        convert_HETATM_to_ATOM,
        # partial(pdb_fixinsert.run, option_list=[]),
        ###
        partial(remove_unsupported_hetatm, user_defined=user_supported_residues),  # noqa: E501
        partial(remove_unsupported_atom, user_defined=user_supported_residues),
        ##
        partial(wdry_pdb_reatom, starting_value=1),
        partial(wdry_pdb_reres, starting_resid=1),
        partial(wdry_pdb_tidy, strict=True),
        ##
        wdry_rstrip,
        ]

    whole_pdb_processing_steps = [
        solve_no_chainID_no_segID,
        ]

    processed_combined_steps = [
        correct_equal_chain_segids,
        ]

    # perform the individual processing steps
    # this list contains a list of lines for each structure
    processed_individually = [
        list(chainf(structure, *line_by_line_processing_steps, dry=dry))
        for structure in structures
        ]

    processed_individually_whole = [
        list(chainf(structure, *whole_pdb_processing_steps))
        for structure in processed_individually
        ]

    processed_combined = chainf(structures, *processed_combined_steps)

    if save_output:
        try:  # in case paths_or_lines is paths
            o_paths = (
                Path(f.parent, f.stem + osuffix).with_suffix(f.suffix)
                for f in map(Path, paths_or_lines)
                )
        except TypeError:  # happens when paths_or_lines was lines
            o_paths = map(
                'processed_{}.pdb'.format,
                range(1, len(structures) + 1),
                )

        for out_p, lines in zip(o_paths, processed_combined):
            out_p.write_text(linesep.join(lines) + linesep)

        return

    else:
        return processed_individually


@_allow_dry("Remove unsupported molecules")
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

    # find a way to report this
    not_allowed_found = set()

    for line in lines:
        if line.startswith(line_startswith):
            residue = line[slc_resname].strip()
            if residue in allowed:
                yield line
            else:
                not_allowed_found.add(residue)
                continue
        else:
            yield line

    return


remove_unsupported_hetatm = partial(
    remove_unsupported_molecules,
    haddock3_defined=supported_HETATM,
    line_startswith='HETATM',
    )

remove_unsupported_atom = partial(
    remove_unsupported_molecules,
    haddock3_defined=supported_ATOM,
    line_startswith='ATOM',
    )


@_allow_dry("Solving chain/seg ID issues.")
def solve_no_chainID_no_segID(lines):
    """
    Solve inconsistencies with chainID and segID.

    If segID is non-existant, copy chainID over segID, and vice-versa.
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


@_allow_dry("Replacing HETATM to ATOM for residue {!r}")
def replace_HETATM_to_ATOM(fhandler, res):
    """."""
    for line in fhandler:
        if line.startswith('HETATM') and line[slc_resname].strip() == res:
            yield 'ATOM  ' + line[6:]
        else:
            yield line


@_allow_dry("Replace residue ATOM/HETATM {!r} to ATOM {!r}")
def replace_residue(fhandler, resin, resout):
    """Replace residue by another and changes HETATM to ATOM if needed."""
    _ = replace_HETATM_to_ATOM(fhandler, res=resin)
    yield from pdb_rplresname.run(_, name_from=resin, name_to=resout)


replace_MSE_to_MET = partial(replace_residue, resin='MSE', resout='MET')
replace_HSD_to_HIS = partial(replace_residue, resin='HSD', resout='HIS')
replace_HSE_to_HIS = partial(replace_residue, resin='HSE', resout='HIS')
replace_HID_to_HIS = partial(replace_residue, resin='HID', resout='HIS')
replace_HIE_to_HIS = partial(replace_residue, resin='HIE', resout='HIS')


@_allow_dry("Add charges to ions.")
def add_charges_to_ions(fhandler):
    """Add charges to ions."""
    for line in fhandler:
        if line.startswith(("ATOM", "ANISOU", "HETATM")):

            # first case
            atom = line[slc_name].strip()
            resname = line[slc_resname].strip()

            # charge is properly defined in resname
            # ignore other fields and write them properly, even if they are
            # already correct
            if resname in supported_ions_resnames:
                new_atom = supported_ions_resnames[resname].atom
                yield line[:12] + new_atom + line[16:76] + new_atom
                continue

            # charge is properly defined in atom name
            # ignore other fields and write them properly, even if they are
            # already correct
            if atom in supported_ions_atoms:
                new_resname = supported_ions_atoms[atom].resname
                yield line[:17] + new_resname + line[20:76] + atom
                continue

            # charge is not defined but atom element is defined
            if atom in supported_ions_elements and atom == resname:
                new_atom = supported_ions_elements[atom].atom
                new_resname = supported_ions_elements[atom].resname
                wmsg = (
                    "Ion {atom!r} automatically set to charge {charge!r}. "
                    "If this is not intended, please edit the PDB manually."
                    )
                log.warning(wmsg)
                yield "".join([
                    line[:12],
                    new_atom,
                    line[16],
                    new_resname,
                    line[20:76],
                    new_atom,
                    ])
                continue

        yield line


def correct_equal_chain_segids(structures):
    """
    Correct chain and segids.

    Parameters
    ----------
    structures : list
    """
    all_chain_ids = set(read_chainids(s) for s in structures)
    remaining_chars = it.cycle(set(_upper_case).difference(all_chain_ids))

    chain_ids = []
    for idx in range(len(structures)):
        chain_id = read_chainids(structures[idx])
        if chain_id in chain_ids:
            _lines = pdb_chain.run(structure, next(remaining_chars))
            structures[idx].clear()
            structures[idx].extend(pdb_chainseg.run(_lines))
        else:
            chain_ids.append(chain_id)

    return  # processed_structures


@_allow_dry("Convert record: {!r} to {!r}.")
def convert_record(fhandler, record, other_record, must_be):
    """Convert ATOM lines to HETATM if needed."""
    for line in fhandler:
        if line.startswith(record):
            resname = line[slc_resname].strip()
            if resname in must_be:
                yield other_record + line[6:]
                continue
        yield line


convert_ATOM_to_HETATM = partial(
    convert_record,
    record="ATOM",
    other_record="HETATM",
    must_be=supported_HETATM,
    )

convert_HETATM_to_ATOM = partial(
    convert_record,
    record="HETATM",
    other_record="ATOM  ",
    must_be=supported_ATOM,
    )
