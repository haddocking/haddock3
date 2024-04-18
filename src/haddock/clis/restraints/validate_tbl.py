"""haddock3-restraints validate_tbl subcommand.

The validate_tbl subcommand validates an input TBL file.

Usage:
    haddock3-restraints validate_tbl <tbl_file> [--silent] [--quick]
"""
from pathlib import Path
from haddock.libs.librestraints import check_parenthesis, validate_tbldata

def add_validate_tbl_arguments(validate_tbl_subcommand):
    """Add arguments to the score subcommand."""
    validate_tbl_subcommand.add_argument(
        "tbl_file",
        type=str,
        help="TBL file to be validated",
        )

    validate_tbl_subcommand.add_argument(
        "--pcs",
        help="PCS mode",
        action='store_true',
        )

    validate_tbl_subcommand.add_argument(
        "--quick",
        help="Check global formatting before going line by line "
        "(opening/closing parenthesis and quotation marks",
        action='store_true',
        )
    
    validate_tbl_subcommand.add_argument(
        "--silent",
        help="Only output errors, do not output TBL file at the end",
        action='store_true',
        )

    return validate_tbl_subcommand


def validate_tbl(tbl_file, pcs, quick=False, silent=False):
    """Get the passive residues."""
    if quick:
        tbldata = open(tbl_file).read()
        # Check the parenthesis and quotation marks opening/closure
        check_parenthesis(tbldata)
    
    if Path(tbl_file).exists():
        tbldata = open(tbl_file).read()
        # Parse and process the restraints
        if silent:
            validate_tbldata(tbldata, pcs)
        else:
            print(validate_tbldata(tbldata, pcs))
    else:
        raise Exception(f"TBL file {tbl_file} does not exist, check the path")
    
    return
