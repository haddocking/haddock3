"""
HADDOCK3 PDB preprocessing client.

Process PDB files for agreement with HADDOCK3 requirements. Follows the
logic implemented in the :py:mod:`haddock.gear.preprocessing`. See
documentation pages for more details.

You can use the `--dry` option to report on the performed changes
without actually performing the changes.

Corrected PDBs are saved to new files named after the `--suffix` option.
Original PDBs are never overwritten, unless `--suffix` is given an empty
string.

You can pass multiple PDB files to the command-line.

Usage::

    haddock-pp file1.pdb file2.pdb
    haddock-pp file1.pdb file2.pdb --suffix _new
    haddock-pp file1.pdb file2.pdb --dry
"""
import argparse
import sys
from pathlib import Path

from haddock.gear.preprocessing import process_pdbs, read_additional_residues
from haddock.libs.libio import add_suffix_to_files, save_lines_to_files
from haddock.libs.libcli import add_output_dir_arg


SUFFIX_DEFAULT = "_processed"

ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    'pdb_files',
    help="Input PDB files.",
    nargs='+',
    )

ap.add_argument(
    '-d',
    '--dry',
    help="Perform a dry run. Informs changes without modifying files.",
    action="store_true",
    )

ap.add_argument(
    '-t',
    '--topfile',
    help="Additional .top files.",
    nargs="*",
    )

ap.add_argument(
    '-s',
    '--suffix',
    help=f"Suffix to output files. Defaults to {SUFFIX_DEFAULT!r}",
    default=SUFFIX_DEFAULT,
    )

add_output_dir_arg(ap)


# client helper functions
def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = vars(load_args(ap))
    # I use this `pop` structure to maintain the unpacking argument in
    # the `main` function because I foresee the `main` to be used from
    # other parts of the software(s)
    main(*cmd.pop("pdb_files"), **cmd)


def maincli():
    """Execute main client."""
    cli(ap, main)


def main(
        *pdb_files,
        dry=False,
        output_directory=None,
        suffix=SUFFIX_DEFAULT,
        topfile=None,
        ):
    """
    Process PDB files.

    Parameters
    ----------
    dry : bool
        Whether to perform a dry test only.

    output_directory : str or ``pathlib.Path``
        The directory where to save the output. Defaults to the current
        working directory.

    suffix : str
        The suffix to append to the new files. Will be added before the
        file extension. Original extension will be kept.

    topfile : str or ``pathlib.Path``
        The path to an additional HADDOCK3 topology file.
    """
    output_directory = output_directory or Path.cwd()
    output_directory.mkdir(parents=True, exist_ok=True)

    new_residues = read_additional_residues(topfile) if topfile else None

    processed_pdbs = process_pdbs(
        *pdb_files,
        dry=dry,
        user_supported_residues=new_residues,
        )

    pdb_names = (Path(output_directory, Path(pdb).name) for pdb in pdb_files)
    out_files = add_suffix_to_files(pdb_names, suffix)
    save_lines_to_files(out_files, processed_pdbs)

    return


if __name__ == '__main__':
    sys.exit(maincli())
