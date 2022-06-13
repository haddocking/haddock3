"""Add functionalities for CLIs."""
from argparse import ArgumentTypeError
from functools import partial

from haddock import version
from haddock.libs.libio import file_exists, folder_exists


arg_file_exist = partial(
    file_exists,
    exception=ArgumentTypeError,
    emsg="File {!r} does not exist or is not a file.",
    )

arg_folder_exist = partial(
    folder_exists,
    exception=ArgumentTypeError,
    emsg="Folder {!r} does not exist.",
    )


def add_version_arg(ap):
    """Add version `-v` argument to client."""
    ap.add_argument(
        "-v",
        "--version",
        help="show version",
        action="version",
        version=f'{ap.prog} - {version}',
        )


def add_rundir_arg(ap):
    """Add run directory option."""
    ap.add_argument(
        "run_dir",
        help="The run directory.",
        type=arg_folder_exist,
        )
