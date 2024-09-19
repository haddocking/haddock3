"""Parse molecular structures in PDB format."""
import os
from functools import partial
from pathlib import Path

from pdbtools.pdb_segxchain import run as place_seg_on_chain
from pdbtools.pdb_splitchain import run as split_chain
from pdbtools.pdb_splitmodel import run as split_model
from pdbtools.pdb_tidy import run as tidy_pdbfile

from haddock.core.supported_molecules import supported_residues
from haddock.core.typing import (
    Callable,
    FilePath,
    FilePathT,
    Iterable,
    Optional,
    Union,
    )
from haddock.libs.libio import working_directory
from haddock.libs.libutil import get_result_or_same_in_list, sort_numbered_paths


slc_record = slice(0, 6)
slc_serial = slice(6, 11)
slc_name = slice(12, 16)
slc_altloc = slice(16, 17)
slc_resname = slice(17, 20)
slc_chainid = slice(21, 22)
slc_resseq = slice(22, 26)
slc_icode = slice(26, 27)
slc_x = slice(30, 38)
slc_y = slice(38, 46)
slc_z = slice(46, 54)
slc_occ = slice(54, 60)
slc_temp = slice(60, 66)
slc_segid = slice(72, 76)
slc_element = slice(76, 78)
slc_charge = slice(78, 80)


def format_atom_name(atom: str, element: str) -> str:
    """
    Format PDB atom name.

    Further Reading:

    * https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html

    Parameters
    ----------
    atom : str
        The atom name.

    element : str
        The atom element code.

    Returns
    -------
    str
        Formatted atom name.
    """
    # string formats for atom name
    _3 = ' {:<3s}'
    _4 = '{:<4s}'
    # len of element, atom formatting string
    # anything outside this is an error
    _atom_format_dict = {
        # len of element: {len of atom name}
        1: {1: _3, 2: _3, 3: _3, 4: _4},
        2: {1: _4, 2: _4, 3: _4, 4: _4},
        }

    atm = atom.strip()
    len_atm = len(atm)
    len_ele = len(element.strip())

    try:
        return _atom_format_dict[len_ele][len_atm].format(atm)
    except KeyError as err:
        _ = f'Could not format this atom:type -> {atom}:{element}'
        # raising KeyError assures that no context in IDPConfGen
        # will handle it. @joaomcteixeira never handles pure Python
        # exceptions, those are treated as bugs.
        raise KeyError(_) from err


def get_supported_residues(haddock_topology: FilePath) -> list[str]:
    """Read the topology file and identify which data is supported."""
    supported: list[str] = []
    with open(haddock_topology) as input_handler:
        for line in input_handler:
            if "resi" in line[:4].casefold():
                res = line.split()[1]
                supported.append(res)
    return supported


_to_remove = ["REMAR", "CTERB", "CTERA", "NTERA", "NTERB", "CONECT"]

_to_rename = {
    "HSD": "HIS",
    "HSE": "HIS",
    "HID": "HIS",
    "HIE": "HIS",
    "WAT ": "TIP3",
    " 0.00969": " 0.00   ",
    }

_to_keep = list(supported_residues)


def split_ensemble(pdb_file_path: Path,
                   dest: Optional[FilePath] = None) -> list[Path]:
    """
    Split a multimodel PDB file into different structures.

    Parameters
    ----------
    dest : str or pathlib.Path
        Destination folder.
    """
    if dest is None:
        dest = Path.cwd()
    assert pdb_file_path.is_file(), pdb_file_path
    with open(pdb_file_path) as input_handler:
        with working_directory(dest):
            split_model(input_handler)

    return sort_numbered_paths(*get_new_models(pdb_file_path))


def split_by_chain(pdb_file_path: FilePath) -> list[Path]:
    """Split a PDB file into multiple structures for each chain."""
    abs_path = Path(pdb_file_path).resolve().parent.absolute()
    with open(pdb_file_path) as input_handler:
        with working_directory(abs_path):
            split_chain(input_handler)

    return get_new_models(pdb_file_path)


def tidy(pdb_file_path: FilePath, new_pdb_file_path: FilePath) -> None:
    """Tidy PDB structure."""
    abs_path = Path(pdb_file_path).resolve().parent.absolute()
    with open(pdb_file_path) as input_handler:
        with working_directory(abs_path):
            with open(new_pdb_file_path, "w") as output_handler:
                for line in tidy_pdbfile(input_handler):
                    output_handler.write(line)


def swap_segid_chain(pdb_file_path: FilePath,
                     new_pdb_file_path: FilePath) -> None:
    """Add to the Chain ID column the found Segid."""
    abs_path = Path(pdb_file_path).resolve().parent.absolute()
    with open(pdb_file_path) as input_handler:
        with working_directory(abs_path):
            with open(new_pdb_file_path, "w") as output_handler:
                for line in place_seg_on_chain(input_handler):
                    output_handler.write(line)


def sanitize(
        pdb_file_path: FilePathT,
        overwrite: bool = True,
        custom_topology: Optional[FilePath] = None) -> Union[FilePathT, Path]:
    """Sanitize a PDB file."""
    if custom_topology:
        custom_res_to_keep = get_supported_residues(custom_topology)
        _to_keep.extend(custom_res_to_keep)

    good_lines: list[str] = []
    with open(pdb_file_path) as input_handler:
        for line in input_handler:
            line = line.rstrip(os.linesep)
            # Ignoring lines containing any tag from __to_remove
            if not any([tag in line for tag in _to_remove]):
                for tag, new_tag in _to_rename.items():
                    line = line.replace(tag, new_tag)
                # check if this residue is known
                res = line[17:20].strip()
                if res and res in _to_keep:
                    good_lines.append(line)
        if len(good_lines) > 0 and good_lines[-1] != "END":
            good_lines.append("END")

    if overwrite:
        with open(pdb_file_path, "w") as output_handler:
            for line in good_lines:
                output_handler.write(line + os.linesep)
        return pdb_file_path

    basename = Path(pdb_file_path)
    new_pdb_file = Path(f"{basename.stem}_cleaned{basename.suffix}")
    new_pdb_file.write_text(os.linesep.join(good_lines) + os.linesep)
    return new_pdb_file


def identify_chainseg(pdb_file_path: FilePath,
                      sort: bool = True) -> tuple[list[str], list[str]]:
    """Return segID OR chainID."""
    segids: list[str] = []
    chains: list[str] = []
    with open(pdb_file_path) as input_handler:
        for line in input_handler:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    segid = line[72:76].strip()[:1]
                except IndexError:
                    segid = ""
                try:
                    chainid = line[21].strip()
                except IndexError:
                    chainid = ""

                if segid:
                    segids.append(segid)
                if chainid:
                    chains.append(chainid)
                
                if not segid and not chainid:
                    raise ValueError(
                        f"Could not identify chainID or segID in pdb {pdb_file_path}, line {line}"
                        )

    if sort:
        segids = sorted(list(set(segids)))
        chains = sorted(list(set(chains)))
    else:
        segids = list(set(segids))
        chains = list(set(chains))
    return segids, chains


def get_new_models(pdb_file_path: FilePath) -> list[Path]:
    """
    Get new PDB models if they exist.

    If no new models are found, return the original path within a list.
    """
    new_models = get_result_or_same_in_list(
        get_pdb_file_suffix_variations,
        pdb_file_path,
        )
    return new_models


def get_pdb_file_suffix_variations(file_name: FilePath,
                                   sep: str = "_") -> list[Path]:
    """
    List suffix variations of a PDB file in the current path.

    If `file.pdb` is given, and files `file_1.pdb`, `file_2.pdb`, exist
    in the folder, those will be listed.

    Parameters
    ----------
    file_name : str or Path
        The name of the file with extension.

    sep : str
        The separation between the file base name and the suffix.
        Defaults to "_".

    Returns
    -------
    list
        List of Paths with the identified PBD files.
        If no files are found return an empty list.
    """
    basename = Path(file_name)
    return list(Path(".").glob(f"{basename.stem}{sep}*{basename.suffix}"))


def read_RECORD_section(
        lines: Iterable[str],
        section_slice: slice,
        func: Callable[[Iterable[str]], Iterable[str]] = set) -> Iterable[str]:
    """
    Create a set of observations from a section of the ATOM line.

    Returns
    -------
    set
        A set of the observations.
    """
    # read the chain ID
    chainids = func(
        the_line
        for line in lines
        if line.startswith(('ATOM', 'HETATM')) and (the_line := line[section_slice].strip())  # noqa: E501
        )
    return chainids


read_chainids = partial(read_RECORD_section, section_slice=slc_chainid, func=list)  # noqa: E501
read_segids = partial(read_RECORD_section, section_slice=slc_segid, func=list)
