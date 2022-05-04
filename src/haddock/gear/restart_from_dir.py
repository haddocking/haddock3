"""Restart from directory gear."""
import shutil
from pathlib import Path

from haddock.core.defaults import MODULE_IO_FILE
from haddock.gear.zerofill import zero_fill
from haddock.libs.libontology import ModuleIO
from haddock.modules import get_module_steps_folders


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
    parent = Path(folder).resolve()
    previous = get_module_steps_folders(folder)

    io = ModuleIO()
    previous_io = Path(parent, previous[0], MODULE_IO_FILE)
    io.load(previous_io)

    return len(io.output)


def renum_step_folders(folder):
    """
    Renumber the step folders sequentially in the run directory.

    The content of the files is not modified.
    See :py:func:`rename_step_contents`.

    Example
    -------
    The initial structure::

        folder/
            0_topoaa/
            4_flexref/

    Results in::

        folder/
            0_topoaa/
            1_flexref/
    """
    # these come sorted already
    step_folders = get_module_steps_folders(folder)

    # move folders
    for i, sf in enumerate(step_folders):
        mod_name = sf.split("_")[-1]
        new_name = zero_fill.fill(mod_name, i)
        shutil.move(Path(folder, sf), Path(folder, new_name))

    return step_folders


def rename_step_contents(folder, previous_step_folders):
    """Rename the contents of the step to be up to date with the new folder name."""
    step_folders = get_module_steps_folders(folder)
    for psf, sf in zip(previous_step_folders, step_folders):
        sfp = Path(folder, sf).resolve()
        for file_ in sfp.iterdir():
            text = file_.read_text()
            text.replace(psf, sf)
            file_.write_text(text)
