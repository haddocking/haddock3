"""Clean module output."""
import gzip
import shutil
import tarfile
from multiprocessing import Pool
from pathlib import Path

from haddock import log
from haddock.libs.libio import (
    archive_files_ext,
    compress_files_ext,
    glob_folder,
    remove_files_with_ext,
    )


UNPACK_FOLDERS = []


def clean_output(path, ncores=1):
    """Perform operations after the run."""
    log.info(f"Cleaning output for {str(path)!r} using {ncores} cores.")
    # add any formats generated to
    # `unpack_compressed_and_archived_files` so that the
    # uncompressing routines when restarting the run work.
    files_to_archive = ['seed', 'inp', 'out', 'con']
    for fta in files_to_archive:
        found = archive_files_ext(path, fta)
        if found:
            remove_files_with_ext(path, fta)

    files_to_compress = ['pdb', 'psf']
    for ftc in files_to_compress:
        found = compress_files_ext(path, ftc, ncores=ncores)
        if found:
            remove_files_with_ext(path, ftc)


def unpack_compressed_and_archived_files(folders, ncores):
    """Unpack compressed and archived files in a folders."""
    global UNPACK_FOLDERS
    UNPACK_FOLDERS.clear()

    for folder in folders:
        gz_files = glob_folder(folder, '.gz')
        tar_files = glob_folder(folder, '.tar')

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
            with tarfile.open(tar_file) as fin:
                fin.extractall(folder)

            tar_file.unlink()


def _unpack_gz(gz_file):
    out_file = Path(gz_file.parent, gz_file.stem)

    with gzip.open(gz_file, 'rb') as fin, \
            open(out_file, 'wb') as fout:
        shutil.copyfileobj(fin, fout, 2 * 10**8)

    gz_file.unlink()
