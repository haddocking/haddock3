"""Features to allow run restart from a given step."""
from argparse import ArgumentParser, ArgumentTypeError
from functools import partial
from pathlib import Path

from haddock.libs.libutil import non_negative_int, remove_folder
from haddock.modules import get_module_steps_folders


_help_cli = """Restart the run from a given step. Previous folders from
the selected step onward will be deleted."""


_arg_non_neg_int = partial(
    non_negative_int,
    exception=ArgumentTypeError,
    emsg="Minimum value is 0, {!r} given.",
    )


def add_restart_arg(parser: ArgumentParser) -> None:
    """Add `--restart` option to argument parser."""
    parser.add_argument(
        "--restart",
        type=_arg_non_neg_int,
        default=None,
        help=_help_cli,
        )


def remove_folders_after_number(run_dir: Path, num: int) -> None:
    """
    Remove calculation folder after (included) a given number.

    Example
    -------
    If the following step folders exist:

        00_topoaa
        01_rigidbody
        02_mdref
        03_flexref

    and the number `2` is given, folders `02_` and `03_` will be
    deleted.

    Parameters
    ----------
    run_dir : pathlib.Path
        The run directory.

    num : int
        The number of the folder from which to delete calculation step
        folders. `num` must be non-negative integer, or equivalent
        representation.
    """
    num = _arg_non_neg_int(num)
    previous = get_module_steps_folders(run_dir.resolve())
    # Filters step folders based on their indices
    from_num_folders = [
        folder for folder in previous
        if int(folder.split('_')[0]) >= num
        ]
    # Loop over folders to remove
    for torm_folder in from_num_folders:
        remove_folder(Path(run_dir, torm_folder))
    return
