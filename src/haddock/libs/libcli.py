"""Add functionalities for CLIs."""
from argparse import ArgumentTypeError
from functools import partial

from haddock import version
from haddock.libs.libutil import file_exists


arg_file_exist = partial(
    file_exists,
    exception=ArgumentTypeError,
    emsg="File {!r} does not exist or is not a file.",
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
