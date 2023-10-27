# TODO:
# 1. single_ions may contain PO4 - verify
# 2. the single_ions_elements map seems not needed
"""
Process input PDB files to ensure compatibility with HADDOCK3.

This module checks and modifies PDB files for compatibility with
HADDOCK3. There are three types of checks/modifications:

1. Performed to each PDB line-by-line, in a equal fashion of ``pdb-tools``.
   In fact, this step mostly uses the ``pdb-tools`` package.
2. Performed on each PDB as a whole.
3. Performed on all PDBs together.

Main functions
--------------

* :py:func:`process_pdbs`
* :py:func:`read_additional_residues`

Corrections performed on 1)
---------------------------

The following actions are perfomed sequentially over all PDBs:

#. from ``pdb-tools``: ``pdb_keepcoord``
#. from ``pdb-tools``: ``pdb_tidy`` with ``strict=True``
#. from ``pdb-toos``: ``pdb_element``
#. from ``pdb-tools``: ``pdb_selaltloc``
#. from ``pdb-tools``: ``pdb_pdb_occ`` with ``occupancy=1.00``
#. replace ``MSE`` to ``MET``
#. replace ``HSD`` to ``HIS``
#. replace ``HSE`` to ``HIS``
#. replace ``HID`` to ``HIS``
#. replace ``HIE`` to ``HIS``
#. add_charges_to_ions, see :py:func:`add_charges_to_ions`
#. convert ``ATOM`` to ``HETATM`` for those atoms that should be ``HETATM``.
   Considers the additional residues provided by the user.
   See :py:func:`convert_ATOM_to_HETATM`.
#. convert ``HETATM`` to ``ATOM`` for those atoms that should be ``ATOM``,
#. from ``pdb-toos``: ``pdb_fixinsert``, with ``option_list=[]``.
#. remove unsupported ``HETATM``. Considers residues provided by the user.
#. remove unsupported ``ATOM``. Considers residues provided by the user.
#. from ``pdb-tools``: ``pdb_reatom``, start from ``1``.
#. from ``pdb-tools``: ``pdb_tidy`` with ``strict=True``

Corrections performed on 2)
---------------------------

The following actions are performed sequentially for each PDB:

* :py:func:`models_should_have_the_same_labels`
* :py:func:`solve_no_chainID_no_segID`
* :py:func:`homogenize_chains`

Read the documentation of the above functions for details what they do.

Corrections performed on 3)
---------------------------

The following actions are performed to all PDBs together:

* :py:func:`correct_equal_chain_segids`

Read the documentation of the above functions for details what they do.

When it happens
---------------

The PDB processing step is performed by default when reading the input
molecules and copying them to the `data/` folder inside the run
directory. When PDBs are processed, a copy of the original input PDBs is
also stored in the `data/` folder.

To deactivate this initial PDB processing, set ``skip_preprocess = False``
in the general parameters of the configuration file.

Additional information
----------------------

If you are a developer and want to read more about the history of this
preprocessing module, visit:

https://github.com/haddocking/haddock3/projects/16
"""
import io
import itertools as it
import re
import string
from functools import partial, wraps
from os import linesep
from pathlib import Path

from pdbtools import (
    pdb_chain,
    pdb_chainxseg,
    pdb_element,
    pdb_fixinsert,
    pdb_keepcoord,
    pdb_occ,
    pdb_reatom,
    pdb_rplresname,
    pdb_segxchain,
    pdb_selaltloc,
    pdb_shiftres,
    pdb_tidy,
)

from haddock import log
from haddock.core.exceptions import HaddockError
from haddock.core.supported_molecules import (
    supported_ATOM,
    supported_HETATM,
    supported_non_ions,
    supported_single_ions_atoms_map,
    supported_single_ions_resnames_map,
)
from haddock.core.typing import (
    Any,
    Callable,
    Container,
    Generator,
    Iterable,
    LineIterSource,
    Optional,
    Union,
)
from haddock.libs.libfunc import chainf
from haddock.libs.libio import read_lines
from haddock.libs.libpdb import (
    format_atom_name,
    read_chainids,
    read_segids,
    slc_charge,
    slc_element,
    slc_name,
    slc_resname,
)


# defines chain letters for chain and seg IDs
_ascii_letters = list(string.ascii_uppercase + string.ascii_lowercase)
_CHAINS = it.cycle(_ascii_letters)


class ModelsDifferError(HaddockError):
    """MODELS of the PDB differ in atom labels."""

    pass


def _report(log_msg: str) -> Callable[..., Any]:
    """
    Add report functionality to the function (decorator).

    Functions decorated with `_report` log the difference between the
    input and the output. Decorated functions gain an additional boolean
    parameter `report` to activate or deactivate the report
    functionality; defaults to ``False``.

    Note that a generator decorated with ``_report`` no longer behaves
    as a generator if ``report=True`` is given. Instead, it returns a
    list from the exhausted generator.

    **Important:** Do NOT use ``_report`` with infinite generators,
    such as ``itertools.cycle``.
    """

    def decorator(function: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(function)
        def wrapper(
            lines: Iterable[Any], *args: Any, report: bool = False, **kwargs: Any
        ) -> Any:
            if report:
                in_lines = list(lines)
                result = list(function(in_lines, *args, **kwargs))

                # Here we could use sets to increase speed, but for the size
                # of the systems, we can actually use lists and get a sorted
                # result by default. I tried using difflib from STD but it is
                # just too slow.

                # _ is line
                additions = [_ for _ in result if _ not in in_lines]
                deletions = [_ for _ in in_lines if _ not in result]

                la = len(additions)
                ld = len(deletions)

                add_lines = linesep.join(f"+ {_}" for _ in additions)
                del_lines = linesep.join(f"- {_}" for _ in deletions)

                log_msg_ = log_msg.format(*args, *kwargs.values())
                extended_log = (
                    f"[{log_msg_}] + {la} - {ld} lines",
                    add_lines,
                    del_lines,
                )

                log.info(linesep.join(extended_log))
                return result

            # If report=False, maintain the original behaviour
            else:
                return function(lines, *args, **kwargs)

        return wrapper

    return decorator


def _open_or_give(inputdata: Iterable[LineIterSource]) -> list[list[str]]:
    """
    Adapt input to the functions.

    Used in py:func:`process_pdbs`.

    Homogenizes input by:

    * removing new line characters at the end of the line
    * removing empty lines

    Parameters
    ----------
    inputdata : list
        A **flat** list where in each index it can contain:

        * file objects
        * paths to files
        * strings representing paths
        * lists or tuples of lines

        The above types can be mixed in the input list.

        Files are read to lines in a list. Line separators are stripped.

        Do not provide nested lists with lists containing paths inside
        lists.

    Returns
    -------
    list of list of strings
        Each sublist has the contents of the input in the same order.

    Raises
    ------
    TypeError
        In any other circumstances.
    """

    def get_line(lines: Iterable[str]) -> list[str]:
        """Ignore empty lines."""
        return [line.rstrip(linesep) for line in lines if line]

    lines: list[list[str]] = []
    for idata in inputdata:
        if isinstance(idata, (Path, str)):
            lines.append(get_line(Path(idata).read_text().split(linesep)))

        elif isinstance(idata, io.TextIOBase):
            lines.append(get_line(idata.readlines()))

        elif isinstance(idata, (list, tuple)):
            lines.append(get_line(idata))

        else:
            emsg = f"Unexpected type in `inputdata`: {type(idata)}"
            raise TypeError(emsg)

    return lines


@read_lines
def read_additional_residues(
    lines: Iterable[str], *ignore: Any, **everything: Any
) -> tuple[str, ...]:
    """
    Read additional residues listed in a ``*.top`` filename.

    Expects new residues to be defined as::

        RESIdue XXX
        RESI XXX
        residue XXX

    where, XXX is the new residue name. Does not read ATOM or charge
    information. Reads only the residue name.

    Examples
    --------
    Read directly the file:

    >>> read_additional_residues(fpath)

    Read the lines instead:

    >>> lines = Path(fpath).read_text().split(os.linesep)
    >>> read_additional_residues.original(lines)

    Parameters
    ----------
    fpath : str or pathlib.Path
        The path to the file.

    lines : list of lines
        You can also use this function in the form of
        ``read_additional_residues.original(...)`` and directly give it
        a list containing the lines of the file.

    Returns
    -------
    tuple
        A tuple with the new identified residues names.
    """
    # https://regex101.com/r/1H44kO/1
    res_regex = re.compile(r"^(RESIdue|residue|RESI) ([A-Z0-9]{1,3}).*$")
    residues: list[str] = []
    for line in map(str.strip, lines):
        name = res_regex.findall(line)
        if name:
            residues.append(name[0][1])

    return tuple(residues)


def process_pdbs(
    *inputdata: LineIterSource,
    dry: bool = False,
    user_supported_residues: Optional[Iterable[str]] = None,
) -> list[list[str]]:
    """
    Process PDB file contents for compatibility with HADDOCK3.

    Parameters
    ----------
    inputdata : list of (str, path, list of str [lines], file handler)

        A **flat** list where in each index it can contain:

        * file objects
        * paths to files
        * strings representing paths
        * lists or tuples of lines

        The above types can be mixed in the input list.

        Files are read to lines in a list. Line separators are stripped.

        Do not provide nested lists with lists containing paths inside
        lists.

    dry : bool
        Perform a dry run. That is, does not change anything, and just
        report.

    user_supported_residues : list, tuple, or set
        The new residues that are allowed.

    Returns
    -------
    list of (list of str)
        The corrected (processed) PDB content in the same order as
        ``inputdata``.
    """
    structures = _open_or_give(inputdata)

    # these are the processing or checking functions that should (if needed)
    # modify the input PDB and return the corrected lines.
    # Follows the same style as for pdb-tools.
    # these functions yield line-by-line.
    line_by_line_processing_steps = [
        wrep_pdb_keepcoord,  # also discards ANISOU
        # tidy is important before some other corrections
        wrep_pdb_tidy_strict,
        wrep_pdb_element,
        wrep_pdb_selaltloc,
        partial(wrep_pdb_occ, occupancy=1.00),
        replace_MSE_to_MET,
        replace_HSD_to_HIS,
        replace_HSE_to_HIS,
        replace_HID_to_HIS,
        replace_HIE_to_HIS,
        add_charges_to_ions,
        partial(
            convert_ATOM_to_HETATM,
            residues=set.union(
                supported_HETATM,
                user_supported_residues or set(),
            ),
        ),
        convert_HETATM_to_ATOM,
        partial(wrep_pdb_fixinsert, option_list=[]),
        #####
        partial(
            remove_unsupported_hetatm, user_defined=user_supported_residues
        ),  # noqa: E501
        partial(remove_unsupported_atom),
        ####
        # partial(wrep_pdb_shiftres, shifting_factor=0),
        partial(wrep_pdb_reatom, starting_value=1),
        wrep_pdb_tidy,
        ###
        wrep_rstrip,
    ]

    # these functions take the whole PDB content, evaluate it, and
    # modify it if needed.
    whole_pdb_processing_steps = [
        models_should_have_the_same_labels,
        solve_no_chainID_no_segID,
        homogenize_chains,
    ]

    # these functions take all structures combined, evulate them
    # togehter, and modify them if needed.
    processed_combined_steps = [
        correct_equal_chain_segids,
    ]

    # START THE ACTUAL PROCESSING

    # individual processing (line-by-line)
    result_1 = [
        list(chainf(structure, *line_by_line_processing_steps, report=dry))
        for structure in structures
    ]

    # whole structure processing
    result_2 = [
        list(chainf(structure, *whole_pdb_processing_steps, report=dry))
        for structure in result_1
    ]

    # combined processing
    final_result = chainf(result_2, *processed_combined_steps)

    return final_result


# Functions operating line-by-line

# make pdb-tools reportable
wrep_pdb_chain = _report("pdb_chain")(pdb_chain.run)
wrep_pdb_chainxseg = _report("pbd_segxchain")(pdb_chainxseg.run)
wrep_pdb_element = _report("pdb_element")(pdb_element.run)
wrep_pdb_fixinsert = _report("pdb_fixinsert")(pdb_fixinsert.run)
wrep_pdb_keepcoord = _report("pdb_keepcoord")(pdb_keepcoord.run)
wrep_pdb_occ = _report("pdb_occ")(pdb_occ.run)
wrep_pdb_reatom = _report("pdb_reatom")(pdb_reatom.run)
wrep_pdb_shiftres = _report("pdb_shiftres")(pdb_shiftres.run)
wrep_pdb_rplresname = _report("pdb_rplresname")(pdb_rplresname.run)
wrep_pdb_segxchain = _report("pdb_segxchain")(pdb_segxchain.run)
wrep_pdb_selaltloc = _report("pdb_selaltloc")(pdb_selaltloc.run)
wrep_pdb_tidy = _report("pdb_tidy")(pdb_tidy.run)
wrep_pdb_tidy_strict = _report("pdb_tidy")(partial(pdb_tidy.run, strict=True))
wrep_rstrip = _report("str.rstrip")(
    partial(map, lambda x: x.rstrip(linesep))
)  # noqa: E501


@_report("Replacing HETATM to ATOM for residue {!r}")
def replace_HETATM_to_ATOM(
    fhandler: Iterable[str], res: str
) -> Generator[str, None, None]:
    """
    Replace record `HETATM` to `ATOM` for `res`.

    Do not alter other lines.

    Parameters
    ----------
    fhanlder : file handler or list of lines
        List-like of file lines. Consumes over a ``for`` loop.

    res : str
        Residue name to match for the substitution.

    Yields
    ------
    str
        Yield line-by-line.
    """
    for line in fhandler:
        if line.startswith("HETATM") and line[slc_resname].strip() == res:
            yield "ATOM  " + line[6:]
        else:
            yield line


@_report("Replace residue ATOM/HETATM {!r} to ATOM {!r}")
def replace_residue(
    fhandler: Iterable[str], resin: str, resout: str
) -> Generator[str, None, None]:
    """
    Replace residue by another and changes ``HETATM`` to ``ATOM`` if needed.

    Do not alter other lines.

    Parameters
    ----------
    fhanlder : file handler or list of lines
        List-like of file lines. Consumes over a ``for`` loop.

    resin : str
        Residue name to match for the substitution.

    resout : str
        Name of the new residue. Renames ``resin`` to ``resout``.

    Yields
    ------
    str
        Yield line-by-line.

    See Also
    --------

    * :py:func:`replace_HETATM_to_ATOM`
    * ``pdb_rplresname`` from ``pdb-tools``
    """
    _ = replace_HETATM_to_ATOM(fhandler, res=resin)
    yield from pdb_rplresname.run(_, name_from=resin, name_to=resout)


replace_MSE_to_MET = partial(replace_residue, resin="MSE", resout="MET")
"""
Replace ``MSE`` to ``MET``.

See Also
--------
* :py:func:`replace_residue`
"""

replace_HSD_to_HIS = partial(replace_residue, resin="HSD", resout="HIS")
"""
Replace ``HSD`` to ``HIS``.

See Also
--------
* :py:func:`replace_residue`
"""

replace_HSE_to_HIS = partial(replace_residue, resin="HSE", resout="HIS")
"""
Replace ``HSE`` to ``HIS``.

See Also
--------
* :py:func:`replace_residue`
"""

replace_HID_to_HIS = partial(replace_residue, resin="HID", resout="HIS")
"""
Replace ``HID`` to ``HIS``.

See Also
--------
* :py:func:`replace_residue`
"""

replace_HIE_to_HIS = partial(replace_residue, resin="HIE", resout="HIS")
"""
Replace ``HIE`` to ``HIS``.

See Also
--------
* :py:func:`replace_residue`
"""


@_report("Remove unsupported molecules")
def remove_unsupported_molecules(
    lines: Iterable[str],
    haddock3_defined: Optional[set[str]] = None,
    user_defined: Optional[set[str]] = None,
    line_startswith: Union[str, tuple[str, ...]] = ("ATOM", "HETATM"),
) -> Generator[str, None, None]:
    """
    Remove HADDOCK3 unsupported molecules.

    This function is abstract and you need to provide the set of
    residues supported by HADDOCK3. See parameters.

    Residues not provided in ``haddock3_defined`` and ``user_defined``
    are removed from the PDB lines.

    Other lines are yieled unmodified.

    Parameters
    ----------
    lines : list or list-like
        Lines of the PDB file. This function will consumes lines over a
        ``for`` loop; mind it if you use a generator.

    haddock3_defined : set
        Set of residues supported by HADDOCK3.
        Defaults to ``None``.

    user_defined : set
        An additional set of allowed residues given by the user.
        Defaults to ``None``.

    line_startswith : tuple
        The lines to consider. Defaults to ``("ATOM", "HETATM")``.

    Yields
    ------
    line : str
        Line-by-line.
        Lines for residues not supported are *not* yielded.

    See Also
    --------
    Other functions use this function to create context.

    * :py:func:`remove_unsupported_atom`
    * :py:func:`remove_unsupported_hetatm`
    """
    user_defined = user_defined or set()
    haddock3_defined = haddock3_defined or set()
    allowed = set.union(haddock3_defined, user_defined)

    # find a way to report this
    not_allowed_found: set[str] = set()

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
    line_startswith="HETATM",
)
"""
Remove unsupported molecules in ``HETATM`` lines.

Uses :py:func:`remove_unsupported_molecules` by populating its
``haddock3_define`` and ``line_startswith`` parameters.

See Also
--------
* :py:func:`remove_unsupported_atom`
"""

remove_unsupported_atom = partial(
    remove_unsupported_molecules,
    haddock3_defined=supported_ATOM,
    line_startswith="ATOM",
)
"""
Remove unsupported molecules in ``ATOM`` lines.

Uses :py:func:`remove_unsupported_molecules` by populating its
``haddock3_define`` and ``line_startswith`` parameters.

See Also
--------
* :py:func:`remove_unsupported_hetatm`
"""


@_report("Add charges to ions.")
def add_charges_to_ions(fhandler: Iterable[str]) -> Generator[str, None, None]:
    """
    Add charges to ions according to HADDOCK3 specifications.

    1. Check if charge is correctly defined in residue name.
       If so, yield the line with correct residue name and charge at the
       end.
    2. Check if charge is correctly defined in atom name.
    3. Create charge from element. This might need manual edit in case
       the atom as an unconventional charge.

    Parameters
    ----------
    fhandler : file-hanlder, list, or list-like
        Lines of the PDB file. This function will consumes lines over a
        ``for`` loop; mind it if you use a generator.

    Yields
    ------
    line : str
        Line-by-line: modified ion lines and any other line.
    """
    # list of functions that correct ion entries
    # by order of preference in case a preference is not found.
    # see further
    ion_correction_cases = [
        _process_ion_case_atom,
        _process_ion_case_resname,
        _process_ion_case_element_charge,
    ]

    for line in fhandler:
        if line.startswith(("ATOM", "ANISOU", "HETATM")):
            # get values
            atom = line[slc_name].strip()  # max 4 chars
            resname = line[slc_resname].strip()  # max 3 chars
            element = line[slc_element].strip()  # max 2 chars
            charge = line[slc_charge].strip()  # max 2 chars

            if resname in supported_non_ions:
                yield line
                continue

            # Which of the above fields has information on the charge?
            # If more than one have information on the charge, we give
            # the preference to the atom, resname, and then charge
            func_to_apply = False
            if atom[-1].isdigit():
                func_to_apply = _process_ion_case_atom  # type: ignore
            elif resname[-1].isdigit():
                func_to_apply = _process_ion_case_resname  # type: ignore
            elif element and charge and charge[-1].isdigit():
                func_to_apply = _process_ion_case_element_charge  # type: ignore

            if func_to_apply:
                yield func_to_apply(line)  # type: ignore

            # in case none of the fields has information on the charge,
            # applies the process functions by giving preference to the
            # residue names.
            else:
                for func in ion_correction_cases:
                    try:
                        yield func(line)
                    except Exception:
                        # test the next function
                        continue
                    else:
                        break  # done
                else:
                    yield line  # lines that do not concern to ions

        else:
            yield line  # lines that are not atoms


def _process_ion_case_resname(line: str) -> str:
    """
    Process ion information based on resnames.

    case 1: charge is correctly defined in resname, for example, ZN2.
    In this case, ignore other fields and write ion information from
    scratch even if it's already correct.
    """
    resname = line[slc_resname].strip()  # max 3 chars
    new_atom = supported_single_ions_resnames_map[resname].atoms[0]
    new_element = supported_single_ions_resnames_map[resname].elements[0]
    charge = new_atom[-2:] if len(new_atom) > 2 else "  "

    new_line = (
        line[:12]
        + format_atom_name(new_atom, new_element)
        + line[16]
        + resname.rjust(3, " ")
        + line[20:76]
        + new_element.rjust(2, " ")
        + charge
    )

    return new_line


def _process_ion_case_atom(line: str) -> str:
    """
    Process ion information based on atom names.

    case 2: charge is correctly defined in atom name ignore other fields
    and write them from scratch even if they are already correct.
    """
    element = line[slc_element].strip()  # max 2 chars
    if element == "C":
        # element C can have atom name CA which conflicts carbon alpha
        # and calcium
        raise ValueError("Element is 'C', does not apply to this case.")

    atom = line[slc_name].strip()  # max 4 chars
    new_resname = supported_single_ions_atoms_map[atom].resname
    new_element = supported_single_ions_atoms_map[atom].elements[0]
    charge = atom[-2:] if len(atom) > 2 else "  "

    new_line = (
        line[:12]
        + format_atom_name(atom, new_element)
        + line[16]
        + new_resname.rjust(3, " ")
        + line[20:76]
        + new_element.rjust(2, " ")
        + charge
    )

    return new_line


def _process_ion_case_element_charge(line: str) -> str:
    """
    Process ion information based on element and charge.

    case 3: charge is correctly defined in atom name ignore other fields
    and write them from scratch even if they are already correct.
    """
    element = line[slc_element].strip()  # max 2 chars
    charge = line[slc_charge].strip()  # max 2 chars
    atom = element + charge
    new_resname = supported_single_ions_atoms_map[atom].resname

    new_line = (
        line[:12]
        + format_atom_name(atom, element)
        + line[16]
        + new_resname.rjust(3, " ")
        + line[20:76]
        + element.rjust(2, " ")
        + charge
    )

    return new_line


@_report("Convert record: {!r} to {!r}.")
def convert_record(
    fhandler: Iterable[str], record: str, other_record: str, residues: Container[str]
) -> Generator[str, None, None]:
    """
    Convert on record to another for specified residues.

    For example, replace ``ATOM`` by ``HETATM`` for specific residues.

    Parameters
    ----------
    fhandler : list-like
        Contains lines of file.

    record : str
        The PDB RECORD to match; for example, ``ATOM`` or ``HETATM``.

    other_record : str
        The PDB RECORD to replace with; for example, ``ATOM`` or ``HETATM``.

    residues : list, tuple, or set
        List of residues to replace the record.
    """
    for line in fhandler:
        if line.startswith(record):
            resname = line[slc_resname].strip()
            if resname in residues:
                yield other_record + line[6:]
                continue
        yield line


convert_ATOM_to_HETATM = partial(
    convert_record,
    record="ATOM",
    other_record="HETATM",
    residues=supported_HETATM,
)
"""
Convert ``ATOM`` to ``HETATM`` for HADDOCK3 supported ``HETATM``.

See Also
--------
* :py:data:`haddock.core.supported_molecules.supported_HETATM`
"""

convert_HETATM_to_ATOM = partial(
    convert_record,
    record="HETATM",
    other_record="ATOM  ",
    residues=supported_ATOM,
)
"""
Convert ``HETATM`` to ``ATOM`` for HADDOCK3 supported ``ATOM``.

See Also
--------
* :py:data:`haddock.core.supported_molecules.supported_ATOM`
"""


# Functions operating in the whole PDB


@_report("Solving chain/seg ID issues.")
def solve_no_chainID_no_segID(lines: Iterable[str]) -> Iterable[str]:
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


# TODO: needs to become an option.
# change chain ID and shift the residue of the other chains.
# also needs to be sync with the restraints.
# maybe not an automatic.. maybe used as a CLI
# maybe is a place to have a CHECK instead of a autoprocess
@_report("Homogenizes chains")
def homogenize_chains(lines: list[str]) -> list[str]:
    """
    Homogenize chainIDs within the same PDB.

    If there are multiple chain identifiers in the PDB file, make all
    them equal to the first one.

    ChainIDs are copied to segIDs afterwards.

    Returns
    -------
    list
        The modified lines.
    """
    chainids = read_chainids(lines)
    if len(set(chainids)) > 1:
        return list(
            chainf(
                lines,
                partial(pdb_chain.run, chain_id=chainids[0]),
                pdb_chainxseg.run,
            )
        )
    else:
        return lines


# Functions operating in all PDBs at once


def correct_equal_chain_segids(structures: list[list[str]]) -> list[list[str]]:
    """
    Correct for repeated chainID in the input PDB files.

    Repeated chain IDs are replaced by an upper case character (``[A-Z]``)
    in order.

    Parameters
    ----------
    structures : list of lists of str
        The input data.

    Returns
    -------
    list of lists of str
        The new structures.
    """
    _all_chains = (read_chainids(s) for s in structures)
    # set of all chain IDs present in the input PDBs
    all_chain_ids = set(it.chain.from_iterable(_all_chains))

    # the remaining available chain characters are the A-Z minus the
    # `all_chain_ids`
    remaining_chars = it.cycle(sorted(set(_ascii_letters).difference(all_chain_ids)))

    chain_ids: list[Iterable[str]] = []
    new_structures: list[list[str]] = []
    for lines in structures:
        new_lines: Optional[list[str]] = None

        # read chain IDs from the structure
        chain_id = read_chainids(lines)

        # if chain_id is repeated
        if chain_id in chain_ids:
            new_lines = list(
                chainf(
                    lines,
                    # change the chain ID by a new one
                    partial(pdb_chain.run, chain_id=next(remaining_chars)),
                    # applies chain ID to seg ID as well
                    pdb_chainxseg.run,
                )
            )
        else:
            chain_ids.append(chain_id)

        new_structures.append(new_lines or lines)

    if len(new_structures) != len(structures):
        raise AssertionError("Number of lines differ. This is a bug!")
    return new_structures


# this id bad
@_report("Check models are the same")
def models_should_have_the_same_labels(lines: Iterable[str]) -> Iterable[str]:
    """
    Confirm models have the same labels.

    In an ensemble of structures, where the PDB file has multiple MODELS,
    all models should have the same labels; hence the same number and
    typ of atoms.

    Parameters
    ----------
    lines : list of strings.
        List containing the lines of the PDB file. Must NOT be a generator.

    Returns
    -------
    list
        The original ``lines`` in case no errors are found.

    Raises
    ------
    ModelsDifferError
        In case MODELS differ. Reports on which models differ.
    """
    # searchers for the first MODEL line. If found, break the loop
    # and continue to the rest of the function.
    #
    # if not found, return the same input lines.
    for line in lines:
        if line.startswith("MODEL"):
            break
    else:
        return lines

    # captures all the models
    models: dict[Optional[int], set[str]] = {}
    new_model: list[str] = []
    new_model_id = None
    for line in lines:
        if line.startswith("MODEL"):
            if new_model_id is not None:
                models[new_model_id] = set(new_model)
                new_model.clear()
            new_model_id = int(line[10:14])

        elif line.startswith(("ATOM", "HETATM")):
            new_model.append(line[12:27])
    else:
        models[new_model_id] = set(new_model)
        new_model.clear()

    # check if all MODELS are equal, performing all vs all comparison
    keys = list(models.keys())
    first_key = keys[0]
    for model_num in keys[1:]:
        if models[model_num] != models[first_key]:
            emsg = f"Labels in MODEL {model_num} differ from MODEL {first_key}."
            raise ModelsDifferError(emsg)

    return lines
