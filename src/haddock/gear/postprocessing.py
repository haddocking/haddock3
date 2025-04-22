"""Tools for post-processing haddock3 runs."""

import os
import shutil
import tarfile

from haddock import log
from haddock.clis.cli_analyse import ANA_FOLDER
from haddock.core.typing import Optional


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
        shutil.rmtree(run_dir)

    return run_archive_fname, analysis_archive_fname