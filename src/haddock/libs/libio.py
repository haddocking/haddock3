"""Lib I/O."""
import contextlib
import glob
import gzip
import os
import tarfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import yaml

from haddock import log
from haddock.libs.libontology import PDBFile
from haddock.libs.libutil import sort_numbered_paths


def clean_suffix(ext):
    """
    Remove the preffix dot of an extension if exists.

    Parameters
    ----------
    ext : str
        The extension string.

    Examples
    --------
    >>> clean_suffix('.pdb')
    'pdb'

    >>> clean_suffix('pdb')
    'pdb'
    """
    return ext.lstrip(r'.')


def dot_suffix(ext):
    """
    Add the dot preffix to an extension if missing.

    Parameters
    ----------
    ext : str
        The extension string.

    Examples
    --------
    >>> clean_suffix('.pdb')
    '.pdb'

    >>> clean_suffix('pdb')
    '.pdb'
    """
    return '.' + clean_suffix(ext)


def read_lines(func):
    """
    Open the file and read lines for the decorated function.

    Send to the decorated function the lines of the file in the form
    of list.
    """
    def wrapper(fpath, *args, **kwargs):
        lines = Path(fpath).read_text().split(os.linesep)
        return func(lines, *args, **kwargs)
    # manual wrapping for displaying documentation properly
    wrapper.original = func
    wrapper.__doc__ = func.__doc__
    wrapper.__name__ = func.__name__
    return wrapper


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


def open_files_to_lines(*files):
    """
    Open files to lines.

    New-lines are stripped.

    Returns
    -------
    list of lists of strings
        The lines of the files.
        Input order is maintained.
    """
    f_paths = map(Path, files)
    return [f.read_text().split(os.linesep) for f in f_paths]


def save_lines_to_files(files, lines):
    """
    Save a list of list of lines to files.

    The first list of strings in `lines` will be saved in the first file
    of `files`, and so on.

    Lines are saved using `pathlib.Path.write_text` function.

    Parameters
    ----------
    files : list
        The list of file names to save.

    lines : list of lists of str
        A list containing lists of lines that are the file contents.
        Must be synched with `files`.
    """
    for file_, content in zip(files, lines):
        Path(file_).write_text(os.linesep.join(content) + os.linesep)

    return


def add_suffix_to_files(files, suffix):
    """
    Add a suffix to file paths.

    Yields
    ------
    pathlib.Path objects
        Exhausts when files exhaust.
    """
    for file_ in files:
        p = Path(file_)
        folder = p.parent
        psuffix = p.suffix
        name = p.stem + suffix + psuffix
        path = Path(folder, name)
        yield path


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


def compress_files_ext(path, ext, ncores=1, **kwargs):
    """
    Compress all files with same extension in folder to `.gz`.

    Do not archive the files in TAR, only compress files individually.

    Parameters
    ----------
    path : str or :external:py:class:`pathlib.Path`
        The folder containing the files.

    ext : str
        The extension of the files.

    **kwargs : anything
        Arguments passed to :py:func:`gzip_files`.

    Returns
    -------
    bool
        ``True`` if files with ``ext`` were found and the compressed
        `.gz` files created.

        ``False`` if no files with ``ext`` were found and, hence, the
        `.gz` files were not created.
    """
    files = glob_folder(path, ext)
    gzip_ready = partial(gzip_files, **kwargs)
    if files:
        with Pool(ncores) as pool:
            imap = pool.imap_unordered(gzip_ready, files)
            for _ in imap:
                pass
        return True
    return False


def gzip_files(file_, block_size=None, compresslevel=9, remove_original=False):
    """
    Gzip a file.

    Parameters
    ----------
    file_ : str or :external:py:class:`pathlib.Path`
        The path to the file to compress.

    block_size : int
        The block size to treat per cycle. Defaults to 200MB (2*10**8
        (2*10**8).

    compresslevel : int
        The compress level. Defaults to 9.
    """
    if block_size is None:
        block_size = 2 * 10**8

    gfile = str(file_) + '.gz'
    with \
            open(file_, 'rb') as fin, \
            gzip.open(gfile, mode='wb', compresslevel=compresslevel) as gout:

        content = fin.read(block_size)  # read the first
        while content:
            gout.write(content)
            content = fin.read(block_size)

    if remove_original:
        Path(file_).unlink()


def archive_files_ext(path, ext, compresslevel=9):
    """
    Archive all files with same extension in folder.

    Parameters
    ----------
    path : str or :external:py:class:`pathlib.Path`
        The folder containing the files.

    ext : str
        The extension of the files.

    compresslevel : int
        The compression level.

    Returns
    -------
    bool
        ``True`` if files with ``ext`` were found and the Zip files created.
        ``False`` if no files with ``ext`` were found and, hence, the
        Zip files was not created.
    """
    files = glob_folder(path, clean_suffix(ext))

    ext = clean_suffix(ext)

    if files:
        with tarfile.open(
                Path(path, f'{ext}.tgz'),
                mode='w:gz',
                compresslevel=compresslevel,
                ) as tarout:

            for file_ in files:
                tarout.add(file_, arcname=file_.name)

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
    ext = f'*{dot_suffix(ext)}'
    files = glob.glob(str(Path(folder, ext)))
    return sort_numbered_paths(*list(map(Path, files)))


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
    files = sort_numbered_paths(*glob_folder(folder, ext))
    # if there are no files, the for loop  won't run.
    for file_ in files:
        log.debug(f'removing: {file_}')
        file_.unlink()


def folder_exists(
        path,
        exception=ValueError,
        emsg="The folder {!r} does not exist or is not a folder.",
        ):
    """
    Assert if a folder exist.

    Parameters
    ----------
    path : str or pathlib.Path
        The path to the folder.

    exception : Exception
        The Exception to raise in case `path` is not file or does not
        exist.

    emsg : str
        The error message to give to `exception`. May accept formatting
        to pass `path`.

    Returns
    -------
    pathlib.Path
        The Path representation of the input ``path`` if condition is
        true.

    Raises
    ------
    Exception
        Any exception that pathlib.Path can raise.
    """
    p = Path(path)

    valid = [p.exists, p.is_dir]

    if all(f() for f in valid):
        return p

    # don't change to f-strings, .format has a purpose
    raise exception(emsg.format(str(path)))


def file_exists(
        path,
        exception=ValueError,
        emsg="`path` is not a file or does not exist",
        ):
    """
    Assert if file exist.

    Parameters
    ----------
    path : str or pathlib.Path
        The file path.

    exception : Exception
        The Exception to raise in case `path` is not file or does not
        exist.

    emsg : str
        The error message to give to `exception`. May accept formatting
        to pass `path`.

    Returns
    -------
    pathlib.Path
        The Path representation of the input ``path`` if condition is
        true.

    Raises
    ------
    Exception
        Any exception that pathlib.Path can raise.
    """
    p = Path(path)

    valid = [p.exists, p.is_file]

    if all(f() for f in valid):
        return p

    # don't change to f-strings, .format has a purpose
    raise exception(emsg.format(str(path)))
