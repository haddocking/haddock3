"""Add functionalities for CLIs."""
from argparse import Action, ArgumentParser, ArgumentTypeError
from functools import partial
from pathlib import Path

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


def add_version_arg(ap: ArgumentParser) -> None:
    """Add version `-v` argument to client."""
    ap.add_argument(
        "-v",
        "--version",
        help="show version",
        action="version",
        version=f"{ap.prog} - {version}",
    )


def add_rundir_arg(ap: ArgumentParser) -> None:
    """Add run directory option."""
    ap.add_argument(
        "run_dir",
        help="The run directory.",
        type=arg_folder_exist,
    )


def add_ncores_arg(ap: ArgumentParser) -> None:
    """Add number of cores option."""
    ap.add_argument(
        "-n",
        "--ncores",
        dest="ncores",
        help=(
            "The number of threads to use. Uses 1 if not specified. "
            "Uses all available threads if `-n` is given. Else, uses the "
            "number indicated, for example: `-n 4` will use 4 threads."
        ),
        type=int,
        default=1,
        const=None,
        nargs="?",
    )


def add_output_dir_arg(ap: ArgumentParser) -> None:
    """Add output dir argument."""
    ap.add_argument(
        "-odir",
        "--output-directory",
        dest="output_directory",
        help="The directory where to save the output.",
        type=Path,
        default=Path.cwd(),
    )


class _ParamsToDict(Action):
    """
    Convert command-line parameters in an argument to a dictionary.

    Example
    -------

    Where ``-x`` is an optional argument of the command-line client
    interface.

        >>> par1 1 par2 'my name' par3 [1,2,3] par4 True
        >>> {'par1': 1, 'par2': 'my name', 'par3': [1, 2, 3]}

    """

    def __call__(self, parser, namespace, ivalues, option_string=None):
        """Execute."""
        params = ivalues[::2]
        values = ivalues[1::2]

        if len(params) != len(values):
            raise parser.error(
                "The parameters and value pairs " "do not match for argument `-p`"
            )

        param_dict = {}
        for k, v in zip(params, values):
            try:
                param_dict[k] = v
            except (ValueError, TypeError, SyntaxError):
                raise parser.error(f"Parameter {k} with invalid value {v}")

        setattr(namespace, self.dest, param_dict)
