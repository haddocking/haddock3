"""
Clean workflow steps' output.

This module concerns removing unnecessary files, compressing, and
archiving files with the same extension to reduce space and stress when
listing files in the modules' step folders.

The two main functions of this module are:

* :py:func:`clean_output`
* :py:func:`unpack_compressed_and_archived_files`

See also the command-line clients ``haddock3-clean`` and
``haddock3-unpack``.
"""
import gzip
import shutil
import tarfile

from functools import partial
from multiprocessing import Pool
from pathlib import Path

from haddock import log
from haddock.core.typing import FilePath, FilePathT, Iterable
from haddock.libs.libio import (
    archive_files_ext,
    compress_files_ext,
    glob_folder,
    remove_files_with_ext,
    )


UNPACK_FOLDERS: list[FilePath] = []


def clean_output(path: FilePath, ncores: int = 1) -> None:
    """
    Clean the output of step folders.

    This functions performs file archiving and file compressing
    operations. Files with extension ``seed``, ``inp``, ``out``, and
    ``con`` are compressed and archived into ``.tgz`` files. The
    original files are deleted.

    Files with ``.pdb`` and ``.psf`` extension are compressed to `.gz`
    files.

    Parameters
    ----------
    path : str or pathlib.Path
        The path to clean. Should point to a folder from a workflow step.

    ncores : int
        The number of cores.
    """
    log.info(f"Cleaning output for {str(path)!r} using {ncores} cores.")
    # add any formats generated to
    # `unpack_compressed_and_archived_files` so that the
    # uncompressing routines when restarting the run work.

    # Files to delete (all)
    file_to_delete_all = (
        "mpi.pkl",
        )
    for extension in file_to_delete_all:
        for file_ in glob_folder(path, extension):
            Path(file_).unlink()

    # Files to delete
    # deletes all except the first one
    # (keeping one for debugging purposes)
    files_to_delete = (
        ".inp",
        ".inp.gz",
        ".out",
        ".out.gz",
        ".job",
        ".err",
        )

    for extension in files_to_delete:
        flist = glob_folder(path, extension)
        for file_ in flist[1:]:
            Path(file_).unlink()

    # files to archive (all files in single .gz)
    files_to_archive = (
        ".seed",
        ".seed.gz",
        ".con",
        )

    archive_ready = partial(_archive_and_remove_files, path=path)
    _ncores = min(ncores, len(files_to_archive))
    with Pool(_ncores) as pool:
        imap = pool.imap_unordered(archive_ready, files_to_archive)
        for _ in imap:
            pass

    # files to compress in .gz
    files_to_compress = (
        ".inp",
        ".out",
        ".pdb",
        ".psf",
        ".cnserr",
        )

    for ftc in files_to_compress:
        found = compress_files_ext(path, ftc, ncores=ncores)
        if found:
            remove_files_with_ext(path, ftc)


def _archive_and_remove_files(fta: str, path: FilePath) -> None:
    found = archive_files_ext(path, fta)
    if found:
        remove_files_with_ext(path, fta)


# eventually this function can be moved to `libs.libio` in case of future need.
def unpack_compressed_and_archived_files(
        folders: Iterable[FilePathT],
        ncores: int = 1,
        dec_all: bool = False,
        ) -> None:
    """
    Unpack compressed and archived files in a folders.

    Works on `.gz` and `.tgz` files.

    Registers folders in :py:data:`UNPACK_FOLDERS` where compressed
    and archived files were found.

    Parameters
    ----------
    folders : list
        List of folders to operate.

    ncores : int
        The number of cores.
    """
    global UNPACK_FOLDERS
    UNPACK_FOLDERS.clear()

    files_to_decompress = [
        ".pdb.gz",
        ".psf.gz",
        ".seed.gz",
        ]

    for folder in folders:
        gz_files: list[Path] = []
        for file_to_dec in files_to_decompress:
            gz_files.extend(list(glob_folder(folder, file_to_dec)))

        if dec_all:
            for dec_all_files in ('.inp.gz', '.out.gz'):
                gz_files.extend(list(glob_folder(folder, dec_all_files)))

        tar_files = glob_folder(folder, '.tgz')

        if gz_files or tar_files:
            # register the folders that where unpacked
            # this is useful for some functionalities of haddock3
            # namely the `--extend-run` option.
            UNPACK_FOLDERS.append(folder)

        if gz_files:  # avoids creating the Pool if there are no .gz files
            with Pool(ncores) as pool:
                imap = pool.imap_unordered(_unpack_gz, gz_files)
                for _ in imap:
                    pass

        for tar_file in tar_files:
            from haddock.libs.libio import extract_files_flat
            extract_files_flat(tar_file, folder)

            tar_file.unlink()


def _unpack_gz(gz_file: Path) -> None:
    out_file = Path(gz_file.parent, gz_file.stem)

    with gzip.open(gz_file, "rb") as fin, \
            open(out_file, "wb") as fout:
        shutil.copyfileobj(fin, fout, 2 * 10**8)

    gz_file.unlink()


def update_unpacked_names(prev: Iterable[FilePath], new: Iterable[FilePath],
                          original: list[FilePath]) -> None:
    """
    Update the unpacked path names.

    Sometimes the step folders are renamed to ajust their index number.
    Such operation happens after the output data is unpacked. This module,
    :py:mod:`haddock.gear.clean_steps`, keeps registry of the folders
    unpacked to the correct funtioning of the `extend_run` module.

    Given the : list[FilePath] names and the new names of the step folders,
    this function updates them in the storing list.

    Examples
    --------
    >>> : list[FilePath] = ['0_topoaa', '4_flexref']
    >>> prev = ['0_topoaa', '4_flexref', '5_seletopclusts']
    >>> new = ['0_topoaa', '1_flexref', '5_seletopclusts']
    >>> update_unpacked_names(prev, new, original)
    >>> original
    ['0_topoaa', '1_flexref']

    This function only evaluate the name of the last folder. And
    maintains the type in the ``original`` list.

    >>> original = ['0_topoaa', Path('4_flexref'), '5_seletopclusts']
    >>> prev = ['0_topoaa', 'run_dir/4_flexref', '5_seletopclusts']
    >>> new = ['run_dir/0_topoaa', '1_flexref', 'run_dir/2_seletopclusts']
    >>> update_unpacked_names(prev, new, original)
    >>> assert original == ['0_topoaa', Path('1_flexref'), '2_seletopclusts']

    Parameters
    ----------
    prev : list of str or pathlib.Path
        The list of the original names before they were changed.

    new : list of str or pathlib.Path
        The list of the new folder names.

    original : list of pathlib.Path
        The list containing the names to record and which names
        will be changed.

    Returns
    -------
    None
        Edits ``original`` in place.
    """
    prevs = map(Path, prev)
    news = map(Path, new)
    original_names = [Path(o).name for o in original]

    # this is used to keep the original types of values in the ``original`` list

    for prev_, new_ in zip(prevs, news):
        try:
            idx = original_names.index(prev_.name)
        except ValueError:  # not present
            continue
        else:
            _p = original[idx]
            original[idx] = type(_p)(Path(Path(_p).parent, new_.name))
