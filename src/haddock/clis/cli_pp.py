"""HADDOCK3 PDB preprocessing client."""
import argparse
import sys

from haddock.gear.preprocessing import process_pdbs, read_additional_residues
from haddock.libs.libio import add_suffix_to_files, save_lines_to_files


ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument('pdb_files', help="Input PDB files.", nargs='+')
ap.add_argument('-d', '--dry', help="Perform a dry run.", action="store_true")
ap.add_argument('-t', '--topfile', help="Additional .top files.", nargs="*")


# client helper functions
def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(ap, main)


def main(pdb_files, dry=False, topfile=None):
    """Process PDB files."""
    new_residues = read_additional_residues(topfile) if topfile else None

    processed_pdbs = process_pdbs(
        pdb_files,
        dry=dry,
        user_supported_residues=new_residues,
        )

    out_files = add_suffix_to_files(pdb_files)
    save_lines_to_files(processed_pdbs, out_files)

    return


if __name__ == '__main__':
    sys.exit(maincli())
