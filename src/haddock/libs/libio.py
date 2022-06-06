"""Lib I/O."""
import contextlib
import glob
import os
from pathlib import Path
from zipfile import ZipFile

import yaml

from haddock import log
from haddock.libs.libontology import PDBFile


def read_from_yaml(yaml_file):
    """
    Read a YAML file to a dictionary.

    Used internally to read HADDOCK3's default configuration files.

    Parameters
    ----------
    yaml_file : str or Path
        Path to the YAML file.

    Returns
    -------
    dict
        Always returns a dictionary.
        Returns empty dictionary if yaml_file is empty.
    """
    with open(yaml_file, 'r') as fin:
        ycfg = yaml.safe_load(fin)

    # ycfg is None if yaml_file is empty
    # returns an empty dictionary to comply with HADDOCK workflow
    if ycfg is None:
        return {}

    assert isinstance(ycfg, dict), type(ycfg)
    return ycfg


def write_dic_to_file(
        data_dict,
        output_fname,
        info_header="",
        sep="\t",
        ):
    """
    Create a table from a dictionary.

    Parameters
    ----------
    data_dict : dict
        Dictionary to write.
    output_fname : str or Path
        Name of the output file.
    info_header : str
        Header to write before the data.
    """
    header = "\t".join(data_dict.keys())

    if info_header:
        header = info_header + os.linesep + header

    with open(output_fname, "w") as out_fh:
        out_fh.write(header + os.linesep)
        row_l = []
        for element in data_dict:
            value = data_dict[element]
            if isinstance(value, Path):
                row_l.append(str(value))
            elif isinstance(value, PDBFile):
                row_l.append(str(value.rel_path))
            elif isinstance(value, int):
                row_l.append(f"{value}")
            elif isinstance(value, str):
                row_l.append(f"{value}")
            elif value is None:
                row_l.append("-")
            else:
                row_l.append(f"{value:.3f}")
        out_fh.write(sep.join(row_l) + os.linesep)


def write_nested_dic_to_file(
        data_dict,
        output_fname,
        info_header="",
        sep="\t"
        ):
    """
    Create a table from a nested dictionary.

    Parameters
    ----------
    data_dict : dict
        Dictionary to write.
    output_fname : str or Path
        Name of the output file.

    Notes
    -----
    This function is used to write nested dictionaries.
    {int: {key: value}}, the int key will be discarded.
    """
    first_key = list(data_dict.keys())[0]
    header = "\t".join(data_dict[first_key].keys())

    if info_header:
        header = info_header + os.linesep + header

    with open(output_fname, "w") as out_fh:
        out_fh.write(header + os.linesep)
        for row in data_dict:
            row_l = []
            for element in data_dict[row]:
                value = data_dict[row][element]
                if isinstance(value, Path):
                    row_l.append(str(value))
                elif isinstance(value, PDBFile):
                    row_l.append(str(value.rel_path))
                elif isinstance(value, int):
                    row_l.append(f"{value}")
                elif isinstance(value, str):
                    row_l.append(f"{value}")
                elif value is None:
                    row_l.append("-")
                else:
                    row_l.append(f"{value:.3f}")
            out_fh.write(sep.join(row_l) + os.linesep)


# thanks to @brianjimenez
@contextlib.contextmanager
def working_directory(path):
    """Change working directory and returns to previous on exit."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


def compress_file_ext(path, ext):
    """
    Compress all files with same extension in folder.

    Parameters
    ----------
    path : str or :external:py:class:`pathlib.Path`
        The folder containing the files.

    ext : str
        The extension of the files.

    Returns
    -------
    bool
        ``True`` if files with ``ext`` were found and the Zip files created.
        ``False`` if no files with ``ext`` were found and, hence, the
        Zip files was not created.
    """
    files = glob_folder(path, ext)
    if files:
        with ZipFile(Path(path, f'{ext}.zip'), 'w') as zipout:
            for file_ in files:
                zipout.write(file_, arcname=file_.name, compresslevel=9)

        return True

    return False


def glob_folder(folder, ext):
    """
    List files with extention `ext` in `folder`.

    Does NOT perform recursive search.

    Parameters
    ----------
    folder : str
        The path to the folder to investigate.

    ext : str
        The file extention. Can be with or without the dot [.]
        preffix.

    Returns
    -------
    list of Path objects
        SORTED list of matching results.
    """
    ext = f'*{parse_suffix(ext)}'
    files = glob.glob(str(Path(folder, ext)))
    return list(map(Path, files))


def parse_suffix(ext):
    """
    Represent a suffix of a file.

    Examples
    --------
    >>> parse_suffix('.pdf')
    '.pdf'

    >>> parse_suffix('pdf')
    '.pdf'

    Parameters
    ----------
    ext : str
        String to extract the suffix from.

    Returns
    -------
    str
        File extension with leading period.
    """
    return f'.{ext[ext.find(".") + 1:]}'


def remove_files_with_ext(folder, ext):
    """
    Remove files with ``ext`` in folder.

    Parameters
    ----------
    folder : str
        The path to the folder.

    ext : str
        The extention of files to delete. Can be with or without the dot ``.``
        preffix.
    """
    files = glob_folder(folder, ext)
    # if there are no files, the for loop  won't run.
    for file_ in files:
        log.debug(f'removing: {file_}')
        file_.unlink()
