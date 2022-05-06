"""
Copy and start run from copy gear.

Contains the functionalities used in `haddock3-copy` CLI and in
`--start-from-copy` flag for `haddock3` CLI.
"""
from pathlib import Path
from shutil import copytree

from haddock import log
from haddock.core.defaults import MODULE_IO_FILE
from haddock.gear.zerofill import zero_fill
from haddock.libs.libontology import ModuleIO
from haddock.libs.libworkflow import Workflow
from haddock.modules import get_module_steps_folders


START_FROM_COPY_DEFAULT = None


class WorkflowManagerCopy:
    """Read and execute workflows from copy."""

    def __init__(self, workflow_params, start=0, **other_params):
        self.start = start
        self.recipe = Workflow(workflow_params, start=start, **other_params)

    def run(self):
        """High level workflow composer."""
        for step in self.recipe.steps:
            step.execute()


def add_start_from_copy(parser):
    """Add option to start-from-dir."""
    parser.add_argument(
        '--start-from-copy',
        help=(
            "Start a run from a run directory previously prepared with "
            "the `haddock3-copy` CLI. Provide the run directory created "
            "with `haddock3-copy` CLI."
            ),
        default=START_FROM_COPY_DEFAULT,
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


def copy_renum_step_folders(indir, destdir, steps):
    """
    Copy step folders renumbering them sequentially in the run directory.

    The content of the files is not modified.
    See :py:func:`rename_step_contents`.

    py:`gear.zerofill.zero_fill`: must be previously calibrated.

    Example
    -------
    >>> steps = ['0_topoaa', '4_flexref']
    >>> zero_fill.set_zerofill_number(len(steps))
    >>> copy_renum_step_folders('run1', 'newrun', steps)

    The initial structure::

        run1/
            1_rigidbody/
            2_caprieval/
            3_seletop/
            4_flexref/
            (etc...)

    Results in::

        newrun/
            0_topoaa/
            1_flexref/

    Parameters
    ----------
    indir : str or Path
        The input directory where the original steps reside.

    destdir : str or Path
        The output directory where to copy the steps to.

    steps : list of (str or Path)
        The list of the folder names in `indir` to copy.
    """
    step_names = (Path(step).name for step in steps)
    for i, step in enumerate(step_names):
        ori = Path(indir, step)
        _modname = step.split("_")[-1]
        dest = Path(destdir, zero_fill.fill(_modname, i))
        copytree(ori, dest)
        log.info(f"Copied {str(ori)} -> {str(dest)}")


def update_contents_of_new_steps(selected_steps, olddir, newdir):
    """
    Find-replace run directory and step name in step folders.

    Parameters
    ----------
    selected_steps : list of str
        The names of the original steps.

    olddir : str or Path
        The original run directory to be replaced by `newdir`.

    newdir : str or Path
        The new run directory.

    Returns
    -------
    None
        Save files in place.
    """
    olddir = Path(olddir)
    newdir = Path(newdir)
    new_steps = get_module_steps_folders(newdir)
    for psf, ns in zip(selected_steps, new_steps):
        new_step = Path(newdir, ns)
        for file_ in new_step.iterdir():
            text = file_.read_text()
            new_text = text.replace(psf, ns)
            new_text = new_text.replace(olddir.name, newdir.name)
            file_.write_text(new_text)

    log.info("File references updated correctly.")
