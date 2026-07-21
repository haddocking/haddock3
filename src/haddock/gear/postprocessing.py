"""Tools for post-processing haddock3 runs."""

import os
import shutil
import tarfile
import time

from haddock import log
from haddock.clis.cli_analyse import ANA_FOLDER
from haddock.core.typing import Optional


def _robust_rmtree(path: str, retries: int = 3, delay: float = 0.5) -> None:
    """Delete a directory tree, tolerant of non-local filesystem quirks.

    On non-local filesystems (NFS, FUSE object-store mounts such as gcsfuse or
    s3fs, overlayfs in containers, CIFS/SMB) deleting a still-open file does not
    remove it immediately, and directory listings can be eventually consistent.
    As a result ``shutil.rmtree`` may fail with
    ``OSError: [Errno 39] Directory not empty`` even when nothing should remain.

    This retries the removal a few times to absorb transient inconsistency. If
    the tree still cannot be removed, a warning is logged and the leftover
    directory is left in place rather than crashing the run — the ``.tgz``
    archive has already been written, so no data is lost.

    The primary safeguard against this situation is closing any open file
    handles (e.g. the run ``log`` file) *before* calling this function; the
    retries here are a defensive fallback for whatever filesystem the run
    happens to live on.
    """
    for attempt in range(1, retries + 1):
        try:
            shutil.rmtree(path)
            return
        except OSError as err:
            if attempt == retries:
                log.warning(
                    f"Could not fully remove '{path}' after archiving "
                    f"({err}). The archive is intact; leaving the directory "
                    "in place. This is typically caused by a non-local "
                    "filesystem (NFS, gcsfuse/s3fs, overlayfs) holding a file "
                    "open or lagging on deletion."
                )
                return
            time.sleep(delay)


def archive_run(run_dir: str, delete: bool = True) -> tuple[str, Optional[str]]:
    """Create an archive of the haddock3 run directory and analysis.

    Parameters
    ----------
    run_dir : str
        Path to the run directory
    delete : bool, optional
        Should the un-archived directory be deleted?, by default False

    Returns
    -------
    tuple[str, Optional[str]]
        run_archive_fname : str
            Path to the run archive
        analysis_archive_fname : Optional[str]
            Path to the run analysis archive
    """
    log.info("Creating an archive of the run")
    # Start by archiving the run_directory
    run_archive_fname = f"{run_dir}.tgz"
    with tarfile.open(run_archive_fname, "w:gz") as tar:
        tar.add(run_dir, arcname=os.path.basename(run_dir))

    # Archive the analysis directory
    analysis_archive_fname = None
    if os.path.exists(f"{run_dir}/{ANA_FOLDER}"):
        analysis_archive_fname = f"{run_dir}_{ANA_FOLDER}.tgz"
        with tarfile.open(analysis_archive_fname, "w:gz") as tar:
            tar.add(f"{run_dir}/{ANA_FOLDER}", arcname=f"{run_dir}_{ANA_FOLDER}")

    if delete:
        _robust_rmtree(run_dir)

    return run_archive_fname, analysis_archive_fname