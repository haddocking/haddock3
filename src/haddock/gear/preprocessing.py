"""
Processes input PDB files to ensure compatibility with HADDOCK3.


"""
from pathlib import Path

from haddock.libs.libfunc import chainf, reduce_helper


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

    f_paths = map(Path, files)
    o_paths = [
        Path(f.parent, f.stem + osuffix, f.suffix)
        for f in map(Path, files)
        ]

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
    individual_processing_steps = [
        ]

    # these are the checks that are performed considering all the PDBs
    # combined
    combined_checking_steps = [
        ]

    # perform the individual processing steps
    processed_individually = [
        chainf(structure, *individual_processing_steps)
        for structure in structures
        ]

    # perform the combined checking steps
    for check_func in combined_checking_steps:
        check_func(processed_individually)

    return processed_individually
