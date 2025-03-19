"""Lib I/O."""
import contextlib
import glob
import gzip
import os
import stat
import tarfile
import re
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import yaml

from haddock import log
from haddock.core.typing import (
    Any,
    Callable,
    FilePath,
    Generator,
    Iterable,
    Mapping,
    Optional,
    Union,
    )
from haddock.libs.libontology import PDBFile
from haddock.libs.libutil import sort_numbered_paths


def clean_suffix(ext: str) -> str:
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
    return ext.lstrip(r".")


def dot_suffix(ext: str) -> str:
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
    return "." + clean_suffix(ext)


def read_lines(func: Callable[..., Any]) -> Callable[..., Any]:
    """
    Open the file and read lines for the decorated function.

    Send to the decorated function the lines of the file in the form
    of list.
    """

    def wrapper(fpath: FilePath, *args: Any, **kwargs: Any) -> Any:
        lines = Path(fpath).read_text().split(os.linesep)
        return func(lines, *args, **kwargs)

    # manual wrapping for displaying documentation properly
    wrapper.original = func  # type: ignore
    wrapper.__doc__ = func.__doc__
    wrapper.__name__ = func.__name__
    return wrapper


def read_from_yaml(yaml_file: FilePath) -> dict[Any, Any]:
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
    # Check that this yaml file do not contain duplicated parameters
    check_yaml_duplicated_parameters(yaml_file)

    # Load yaml file using the yaml lib
    with open(yaml_file, "r") as fin:
        ycfg = yaml.safe_load(fin)

    # ycfg is None if yaml_file is empty
    # returns an empty dictionary to comply with HADDOCK workflow
    if ycfg is None:
        return {}

    assert isinstance(ycfg, dict), type(ycfg)
    return ycfg


def check_yaml_duplicated_parameters(yaml_fpath: str) -> None:
    """Make sure the provided yaml file do not contain duplicated parameters.

    Parameters
    ----------
    yaml_fpath : str
        Path to a yaml file
    """
    # Build regular expression
    # Note: Understand behavior here -> https://regex101.com/r/AaFHp4/1
    yaml_param_regex = re.compile("^(([A-Za-z0-9]_?)+):")
    # Read content as string
    with open(yaml_fpath, 'r') as filin:
        yaml_content = filin.readlines()
    parsed_param_names: dict[str, int] = {}
    # Loop over lines
    for i, line in enumerate(yaml_content, start=1):
        # Check if new parameter
        if (match := yaml_param_regex.search(line)):
            # Point parameter name
            param_name = match.group(1)
            # Make sure this parameter has not yet been used
            assert param_name not in parsed_param_names.keys(), f"Parameter '{param_name}' in {yaml_fpath} has duplicates: l.{parsed_param_names[param_name]} and l.{i}"  # noqa : E501
            # Hold line were this parameter is, in case of duplication, to help
            parsed_param_names[param_name] = i


def open_files_to_lines(*files: FilePath) -> list[list[str]]:
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


def save_lines_to_files(
        files: Iterable[FilePath],
        lines: Iterable[Iterable[str]]
        ) -> None:
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


def add_suffix_to_files(
        files: Iterable[FilePath],
        suffix: str
        ) -> Generator[Path, None, None]:
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
        data_dict: Mapping[Any, Any],
        output_fname: FilePath,
        info_header: str = "",
        sep: str = "\t",
        ) -> None:
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
        row_l: list[str] = []
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
        data_dict: Mapping[Any, Any],
        output_fname: FilePath,
        info_header: str = "",
        sep: str = "\t",
        ) -> None:
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
            row_l: list[str] = []
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
def working_directory(path: FilePath) -> Generator[None, None, None]:
    """Change working directory and returns to previous on exit."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


def compress_files_ext(
        path: FilePath,
        ext: str,
        ncores: int = 1,
        **kwargs: Any
        ) -> bool:
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


def gzip_files(
        file_: FilePath,
        block_size: Optional[int] = None,
        compresslevel: int = 9,
        remove_original: bool = False,
        ) -> None:
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

    gfile = str(file_) + ".gz"
    with open(file_, "rb") as fin, gzip.open(
            gfile, mode="wb", compresslevel=compresslevel
            ) as gout:
        content = fin.read(block_size)  # read the first
        while content:
            gout.write(content)
            content = fin.read(block_size)

    if remove_original:
        Path(file_).unlink()


def archive_files_ext(path: FilePath, ext: str, compresslevel: int = 9) -> bool:
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
                Path(path, f"{ext}.tgz"),
                mode="w:gz",
                compresslevel=compresslevel,
                ) as tarout:
            for file_ in files:
                tarout.add(file_, arcname=file_.name)

        return True
    return False


def glob_folder(folder: FilePath, ext: str) -> list[Path]:
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
    ext = f"*{dot_suffix(ext)}"
    files = glob.glob(str(Path(folder, ext)))
    return sort_numbered_paths(*(Path(file) for file in files))


def remove_files_with_ext(folder: FilePath, ext: str) -> None:
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
        log.debug(f"removing: {file_}")
        file_.unlink()


def folder_exists(
        path: FilePath,
        exception: type[Exception] = ValueError,
        emsg: str = "The folder {!r} does not exist or is not a folder.",
        ) -> Path:
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
        path: FilePath,
        exception: type[Exception] = ValueError,
        emsg: str = "`path` is not a file or does not exist",
        ) -> Path:
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


def pdb_path_exists(pdb_path: Path) -> tuple[bool, Optional[str]]:
    """
    Check if a pdb path exists.

    If not, checks for the existence of a gzipped pdb file and informs the user
    that the file is gzipped

    Parameters
    ----------
    pdb_path : pathlib.Path
        path to the pdb

    Returns
    -------
    exists : bool
        True if file exists
    msg : str or None
        the error message
    """
    exists, msg = True, None
    if not pdb_path.exists():
        msg = f"PDB file {pdb_path} not found."
        gz_pdb_path = pdb_path.with_suffix(pdb_path.suffix + ".gz")
        if gz_pdb_path.exists():
            msg += f" A compressed file ({gz_pdb_path}) exists though."
            msg += "Use haddock3-unpack to unpack the run."
        exists = False
    return exists, msg


def get_perm(fname: FilePath) -> int:
    """Get permissions of file."""
    # https://stackoverflow.com/questions/6874970
    return stat.S_IMODE(os.lstat(fname)[stat.ST_MODE])


def make_writeable_recursive(path: FilePath) -> None:
    """
    Add writing to a folder, its subfolders and files.

    Parameters
    ----------
    path : str or Path
        The path to add writing permissions.
    """
    # https://stackoverflow.com/questions/6874970
    for root, dirs, files in os.walk(path, topdown=False):
        for dir_ in (os.path.join(root, d) for d in dirs):
            os.chmod(dir_, get_perm(dir_) | stat.S_IWUSR)

        for file_ in (os.path.join(root, f) for f in files):
            os.chmod(file_, get_perm(file_) | stat.S_IWUSR)


def extract_files_flat(tar_path: Union[str, FilePath],
                       dest_path: Union[str, FilePath]) -> None:
    """
    Extract files from a tarball to a destination folder.
    
    Parameters
    ----------
    tar_path : str or Path
        The path to the tarball file.

    dest_path : str or Path
        The path to the destination folder.
    """
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            # Extract only files (skip directories)
            if member.isfile():
                # Modify the member name to remove the directory structure
                member.name = os.path.basename(member.name)
                tar.extract(member, dest_path)
