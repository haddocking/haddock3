"""Gear for ``haddock3-copy`` CLI and `--extend-run`` flag."""
import shutil
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_IO_FILE
from haddock.core.exceptions import HaddockTermination
from haddock.core.typing import (
    Any,
    ArgumentParser,
    FilePath,
    FilePathT,
    Iterable,
    ModuleParams,
    Optional,
    )
from haddock.gear.clean_steps import UNPACK_FOLDERS, clean_output
from haddock.gear.zerofill import zero_fill
from haddock.libs.libontology import ModuleIO
from haddock.libs.libtimer import log_time
from haddock.libs.libworkflow import Workflow, WorkflowManager
from haddock.modules import get_module_steps_folders


EXTEND_RUN_DEFAULT = None


class WorkflowManagerExtend(WorkflowManager):
    """Workflow to extend a run."""

    def __init__(self,
                 workflow_params: ModuleParams,
                 start: Optional[int] = 0,
                 **other_params: Any) -> None:
        self.start = start
        self.recipe = Workflow(workflow_params, start=start, **other_params)
        # terminate is used to synchronize the `clean` option with the
        # `exit` module. If the `exit` module is removed in the future,
        # you can also remove and clean the `terminate` part here.
        self._terminated = 0

    def run(self) -> None:
        """High level workflow composer."""
        for i, step in enumerate(self.recipe.steps, start=0):
            try:
                step.execute()
            except HaddockTermination:
                self._terminated = i
                break

    def clean(self) -> None:
        """Clean the step output."""
        # return compression to the original state
        cwd = Path.cwd().name

        # because this WorkflowManagerExtended has no direct access to the
        # ncores parameters, it needs to take it from the steps.
        ncores: int = max(s.config["ncores"] for s in self.recipe.steps)

        for folder in UNPACK_FOLDERS:
            # temporary hack to get the step folder name.
            # ensures we can work under the CLI working directory
            # or inside the run dir.
            folder_ = str(folder).split(cwd)[1][1:]

            log.info(
                f'Compressing original folder: {folder_!r} because it '
                'was originally compressed.'
                )

            with log_time("cleaning output files took"):
                clean_output(folder_, ncores)

        # apply compression to the new modules
        super().clean(terminated=self._terminated)


def add_extend_run(parser: ArgumentParser) -> None:
    """Add option to ``--extend-run``."""
    parser.add_argument(
        '--extend-run',
        help=(
            "Start a run from a run directory previously prepared with "
            "the `haddock3-copy` CLI. Provide the run directory created "
            "with `haddock3-copy` CLI."
            ),
        default=EXTEND_RUN_DEFAULT,
        type=Path,
        )


def read_num_molecules_from_folder(folder: FilePath) -> int:
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


def copy_renum_step_folders(indir: FilePath, destdir: FilePath,
                            steps: list[FilePathT]) -> list[Path]:
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
            0_topoaa/
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

    Returns
    -------
    list
        The new paths created.
    """
    new_steps: list[Path] = []
    step_names = (Path(step).name for step in steps)
    for i, step in enumerate(step_names):
        ori = Path(indir, step)
        _modname = step.split("_")[-1]
        dest = Path(destdir, zero_fill.fill(_modname, i))
        shutil.copytree(ori, dest)
        log.info(f"Copied {str(ori)} -> {str(dest)}")
        new_steps.append(dest)
    return new_steps


def renum_step_folders(folder: FilePath) -> tuple[list[str], list[str]]:
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

    Returns
    -------
    list
        The list of the original step folder names.

    list
        The list of the new step folder names.
    """
    # these come sorted already
    step_folders = get_module_steps_folders(folder)
    new_names: list[str] = []

    # move folders
    for i, sf in enumerate(step_folders):
        mod_name = sf.split("_")[-1]
        new_name = zero_fill.fill(mod_name, i)
        new_names.append(new_name)
        shutil.move(Path(folder, sf), Path(folder, new_name))

    return step_folders, new_names


def update_contents_of_new_steps(selected_steps: Iterable[str],
                                 olddir: FilePath, newdir: FilePath) -> None:
    """
    Find-replace run references in step folders files.

    Find and replaces (updates) all references to step folders and to
    the old run directory in all files of selected step folders.

    Example
    -------
    >>> update_contents_of_new_steps(
        ['0_topoaa', '1_rigidbody'],
        'run1',
        'run2',
        )

    Parameters
    ----------
    selected_steps : list of str
        The names of the original step folder names that to find in the
        new step folders. This function uses
        :py:func:`haddock.modules.get_module_steps_folders` to find
        the new step folders. ``selected_steps`` must be synchronized
        with the new folders; that is, the names in ``selected_steps``
        must be the old names of the new step folders.

    olddir : str or Path
        The original run directory to be replaced by `newdir`.

    newdir : str or Path
        The new run directory.

    Returns
    -------
    None
        Save files in place.

    See Also
    --------
    :py:func:`haddock.gear.prepare_run.update_step_contents_to_step_names`
    """
    olddir = Path(olddir)
    newdir = Path(newdir)
    new_steps = get_module_steps_folders(newdir)
    for ns in new_steps:
        new_step = Path(newdir, ns)
        for file_ in new_step.iterdir():
            try:
                text = file_.read_text()
            except UnicodeDecodeError as err:
                log.warning(f"Failed to read file {file_}. Error is {err}")
                continue
            for s1, s2 in zip(selected_steps, new_steps):
                text = text.replace(s1, s2)
            text = text.replace(olddir.name, newdir.name)
            file_.write_text(text)

    log.info("File references updated correctly.")
