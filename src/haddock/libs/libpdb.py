"""Parse molecular structures in PDB format."""
import os
from pathlib import Path

from pdbtools.pdb_segxchain import place_seg_on_chain
from pdbtools.pdb_splitchain import split_chain
from pdbtools.pdb_splitmodel import split_model
from pdbtools.pdb_tidy import tidy_pdbfile

from haddock.core.cns_paths import topology_file
from haddock.libs.libio import working_directory
from haddock.libs.libutil import get_result_or_same_in_list, sort_numbered_paths


def get_supported_residues(haddock_topology):
    """Read the topology file and identify which data is supported."""
    supported = []
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

_to_keep = get_supported_residues(topology_file)


def split_ensemble(pdb_file_path, dest=None):
    """
    Split a multimodel PDB file into different structures.

    Parameters
    ----------
    dest : str or pathlib.Path
        Destination folder.
    """
    dest = Path.cwd()
    assert pdb_file_path.is_file(), pdb_file_path
    with open(pdb_file_path) as input_handler:
        with working_directory(dest):
            split_model(input_handler)

    return sort_numbered_paths(*get_new_models(pdb_file_path))


def split_by_chain(pdb_file_path):
    """Split a PDB file into multiple structures for each chain."""
    abs_path = Path(pdb_file_path).resolve().parent.absolute()
    with open(pdb_file_path) as input_handler:
        with working_directory(abs_path):
            split_chain(input_handler)

    return get_new_models(pdb_file_path)


def tidy(pdb_file_path, new_pdb_file_path):
    """Tidy PDB structure."""
    abs_path = Path(pdb_file_path).resolve().parent.absolute()
    with open(pdb_file_path) as input_handler:
        with working_directory(abs_path):
            with open(new_pdb_file_path, "w") as output_handler:
                for line in tidy_pdbfile(input_handler):
                    output_handler.write(line)


def swap_segid_chain(pdb_file_path, new_pdb_file_path):
    """Add to the Chain ID column the found Segid."""
    abs_path = Path(pdb_file_path).resolve().parent.absolute()
    with open(pdb_file_path) as input_handler:
        with working_directory(abs_path):
            with open(new_pdb_file_path, "w") as output_handler:
                for line in place_seg_on_chain(input_handler):
                    output_handler.write(line)


def sanitize(pdb_file_path, overwrite=True, custom_topology=False):
    """Sanitize a PDB file."""
    if custom_topology:
        custom_res_to_keep = get_supported_residues(custom_topology)
        _to_keep.extend(custom_res_to_keep)

    good_lines = []
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


def identify_chainseg(pdb_file_path, sort=True):
    """Return segID OR chainID."""
    segids = []
    chains = []
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

    if sort:
        segids = sorted(list(set(segids)))
        chains = sorted(list(set(chains)))
    else:
        segids = list(set(segids))
        chains = list(set(chains))
    return segids, chains


def get_new_models(pdb_file_path):
    """
    Get new PDB models if they exist.

    If no new models are found, return the original path within a list.
    """
    new_models = get_result_or_same_in_list(
        get_pdb_file_suffix_variations,
        pdb_file_path,
        )
    return new_models


def get_pdb_file_suffix_variations(file_name, path=None, sep="_"):
    """
    List suffix variations of a PDB file.

    If `file.pdb` is given, and files `file_1.pdb`, `file_2.pdb`, exist
    in the folder, those will be listed.

    Parameters
    ----------
    file_name : str or Path
        The name of the file with extension.

    path : str or pathlib.Path
        Path pointing to a directory where to perform the search.

    sep : str
        The separation between the file base name and the suffix.
        Defaults to "_".

    Returns
    -------
    list
        List of Paths with the identified PBD files.
        If no files are found return an empty list.
    """
    folder = path or Path.cwd()

    if not folder.is_dir():
        raise ValueError(f'{str(folder)!r} should be a directory.')

    basename = Path(file_name)
    return list(folder.glob(f"{basename.stem}{sep}*{basename.suffix}"))
