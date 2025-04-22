"""Parse molecular structures in PDB format."""
import os
from functools import partial
from pathlib import Path

from pdbtools.pdb_segxchain import run as place_seg_on_chain
from pdbtools.pdb_splitchain import run as split_chain
from pdbtools.pdb_splitmodel import run as split_model
from pdbtools.pdb_tidy import run as tidy_pdbfile

from haddock.core.exceptions import SetupError
from haddock.core.supported_molecules import supported_residues
from haddock.core.typing import (
    Callable,
    FilePath,
    FilePathT,
    Iterable,
    Optional,
    TypeVar,
    Union,
    )
from haddock.libs.libio import working_directory, PDBFile
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

# Get list of supported residues to keep when pre-processing
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
        # Terminate file with an END statement
        if len(good_lines) > 0 and good_lines[-1] != "END":
            good_lines.append("END")

    # Check if anything has been kept from this file
    if len(good_lines) == 0:
        raise SetupError(
            f"No coordinates kept after sanitzing {pdb_file_path.name}."
        )

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


def add_TER_on_chain_breaks(
        input_pdb: FilePath,
        output_pdb: FilePath,
        ) -> None:
    """Detect chain breaks and add TER statements between them.

    Parameters
    ----------
    input_pdb : FilePath
        Input PDB filepath with potential chain breaks.
    output_pdb : FilePath
        Output PDB filepath with added TER statements between chain breaks.
    """
    Residue = dict[str, Union[list[str], list[float], str]]
    residueT = TypeVar('residueT', bound=Residue)

    def euclidean_dist(atm1: list[float], atm2: list[float]) -> float:
        """Compute Euclidean distances between two points.

        Parameters
        ----------
        atm1 : list[float]
            Atom 1 coordinates.
        atm2 : list[float]
            Atom 2 coordinates.

        Returns
        -------
        float
            Distance between the two atoms.
        """
        return sum([(c1 - c2) ** 2 for c1, c2 in zip(atm1, atm2)]) ** 0.5

    def detected_chain_break(residue1: residueT, residue2: residueT) -> bool:
        """Detect chain break between two consecutive residues.

        Currently limitted to DNA and protein chains.

        Parameters
        ----------
        residue1 : residueT
            Previous residue (N-1)
        residue2 : residueT
            Current residue (N)

        Returns
        -------
        bool
            True if chain break detected, else False
        """
        # Expected backbone distances set
        backbone_dists = {
            "protein": 3.5,  # very loose !
            "DNA": 4.5,
            }
        # Detect type of residues
        # DNA case
        try:
            # Extract coordinates
            atm1 = residue1["O3'"]
            atm2 = residue2["O5'"]
        except KeyError:
            # Protein case
            try:
                # Extract coordinates
                atm1 = residue1["C"]  # C of -1 residue
                atm2 = residue2["N"]  # N of 0 residue
            # Error in detection of any of the atoms
            except KeyError:
                # Must be a chain break
                # FIXME : This assumes that it is NOT a DNA nor a protein.
                # In the future, if future there is, it may cause an issue if
                # lipids, glycans, other type of entities are added,
                # as this will trigger a TER statement between them.
                return True
            else:
                entity_type = "protein"
        else:
            entity_type = "DNA"
        # Point distance
        upper_dist = backbone_dists[entity_type]
        # Check if distance is within acceptable peptide bond limit
        return euclidean_dist(atm1, atm2) > upper_dist

    def write_residue(fhandler, residue_lines: list[str]) -> None:
        """Writes residues line to file.
        
        Parameters
        ----------
        fhandler : _type_
            File object on which to write the residue lines.
        residue_lines : list[str]
            Residue to write.
        """
        for _ in residue_lines:
            fhandler.write(_)
    
    def write_previous_residue(
            fhandler,
            previous: residueT,
            current: residueT,
            ) -> None:
        """Write residue lines to file, possibly ending by TER.

        Parameters
        ----------
        fhandler : _type_
            File object on which to write the residue lines.
        previous : residueT
            Previous residue (N-1)
        current : residueT
            Current residue (N)
        """
        # Write previous residue
        write_residue(fhandler, previous["lines"])
        # Check if bond observed between -2 and -1
        chain_break = detected_chain_break(
            previous,
            current,
            )
        # If chain break detected or new chain in file
        if chain_break or previous["chain"] != current["chain"]:
            fhandler.write(f"TER{os.linesep}")

    # Initiate parsing variables
    BB_atomnames: tuple[str, str, str, str] = (
        "C", "N",  # for peptide bonds
        "O3'", "O5'",  # for DNA/RNA
        )
    current_resid: tuple[str, str] = ("-", "-", )
    previous_residue: residueT = {"lines": []}
    current_residue: residueT = {"lines": []}
    # Read input file
    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        for _ in fin:
            if _.startswith(("ATOM", "HETATM", )):
                # Extract specific data from coordinates record
                chainid = _[slc_chainid].strip()
                resid = _[slc_resseq].strip()
                atname = _[slc_name].strip()

                # Case when new residue
                if current_resid != (chainid, resid):
                    # Make sure it is not first residue
                    if len(previous_residue["lines"]) > 0:
                        # Write previous residue
                        write_previous_residue(
                            fout,
                            previous_residue,
                            current_residue
                            )
                    # Reset previous residue to current residue
                    previous_residue = current_residue
                    # Initialize new current residue
                    current_resid = (chainid, resid, )
                    current_residue = {"lines": [], "chain": chainid}

                # Hold residue line
                current_residue["lines"].append(_)

                # If atom name is in the backbone atoms
                if atname in BB_atomnames:
                    # Hold backbone atom coordinates in current residue
                    current_residue[atname] = [
                        float(_[slc_x]),
                        float(_[slc_y]),
                        float(_[slc_z]),
                        ]
        # Last previous residue
        write_previous_residue(
            fout,
            previous_residue,
            current_residue
            )
        # Last residue
        write_residue(fout, current_residue["lines"])
        # Write final TER statement
        fout.write(f"TER{os.linesep}")
        # Write END statement
        fout.write(f"END{os.linesep}")


def check_combination_chains(combination: list[PDBFile]) -> list[str]:
    """Check if chain IDs are unique for each pdb in combination."""
    chainid_list: list[str] = []
    for pdb in combination:
        segids, chains = identify_chainseg(pdb.rel_path, sort=False)
        chainsegs = sorted(list(set(segids) | set(chains)))
        # check if any of chainsegs is already in chainid_list
        if any(chainseg in chainid_list for chainseg in chainsegs):
            raise ValueError(
                f"Chain/seg IDs are not unique for pdbs {combination}."
            )
        chainid_list.extend(chainsegs)
    return chainid_list
