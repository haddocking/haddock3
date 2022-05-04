"""Restart from directory gear."""
from pathlib import Path

from haddock.core.defaults import MODULE_IO_FILE, modules_folder_prefix
from haddock.libs.libontology import ModuleIO
from haddock.libs.libutil import glob_step_folders


RESTART_FROM_DIR_DEFAULT = None


def add_restart_from_dir(parser):
    """Add option to restart-from-dir."""
    parser.add_argument(
        '--restart-from-dir',
        help="The run directory to restart from.",
        default=RESTART_FROM_DIR_DEFAULT,
        type=Path,
        )


def read_num_molecules_from_folder(folder):
    """
    Read the number of molecules from the first step folder.

    1. Find the lower indexed step folder in `folder`.
    2. Read the `io.json` file.
    3. Count the number of items in the "output" key.
    4. The above is the number of molecules.

    Parameters
    ----------
    folder : Path or string
        The run directory.

    Returns
    -------
    int
        The number of molecules found.
    """
    previous = glob_step_folders(folder)

    io = ModuleIO()
    previous_io = Path(previous[0], MODULE_IO_FILE)
    io.load(previous_io)

    return len(io.output)
