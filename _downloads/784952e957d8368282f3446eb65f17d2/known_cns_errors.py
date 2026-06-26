"""Detect known/common CNS errors.

Inspired from:
https://github.com/haddocking/haddock25/blob/main/tools/check-error-messages.sh
"""

import gzip
from io import BufferedReader
from pathlib import Path

from haddock.core.exceptions import KnownCNSError
from haddock.core.typing import FilePath, Optional, Union

# Dictionary of known errors
# as key:    How to catch it in the cns.cnserr
# as value:  Message to user
KNOWN_ERRORS = {
    "CHAIN LENGTH FOR SYMMETRY RESTRAINTS DOES NOT MATCH": (
        "Mismatch between chain length for symmetry restraints. "
        "Check your input molecules and symmetry restraints."
        ),
    "NCS-restraints error encountered: Improperly defined non-crystallographic symmetry": (  # noqa : E501
        "Improperly defined non-crystallographic symmetry (NCS). "
        "Check your symmetry restraints definition."
        ),
    "error in SYMMETRY potential, check NOE table": (
        "Check your symmetry restraints definition."
        ),
    "exceeded allocation for NOE-restraints": (
        "Too many distance restraints defined. "
        "Try to reduce this number by checking your definition of active "
        "and passive residues. "
        "Make sure to filter those for solvent accessibility. "
        "Or alternatively increase the nres parameter in the noe statements"
        " in the relevant CNS scripts."
        ),
    "SELRPN error encountered: parsing error": (
        "Check your restraint files."
        ),
    "PARSER error encountered: Encountered too many parsing errors": (
        "Encountered too many parsing errors. "
        "Check your input molecules and symmetry restraints."
        ),
    "XMREAD error encountered:  sectioning of map incompatible with resolution": (  # noqa : E501
        "Check your EM map resolution and sectioning."
        ),
    "ALLHP error encountered: not enough memory available": (
        "Too many distance restraints defined. "
        "Try to reduce this number by checking your definition of active and "
        "passive residues. "
        "Make sure to filter those for solvent accessibility. "
        "Try to decrease the size of your system where possible."
        ),
    "error encountered: missing SCATter definition for SELEcted atoms": (
        "Unsupported atoms/molecules for cryo-EM restraints."
        ),
    "ROTMAT error encountered: rotation vector has zero length": (
        "Check your input parameters and restraints. "
        "Possibly try turning off the sampling of 180 degrees rotation."
        )
    }


def find_cns_errors(cns_out_fpath: FilePath) -> Optional[KnownCNSError]:
    """Detect if a known CNS error is in a cns.cnserr file.

    Parameters
    ----------
    cns_out_fpath : FilePath -> Union[str, Path]
        Path to the cns.cnserr file to check.

    Returns
    -------
    Optional[KnownCNSError]
        An exception for known CNS errors, with its hint on how to solve it!
    """
    # Check for file extension to open it the appropriate way
    if Path(cns_out_fpath).suffix == ".gz":
        file_handle = gzip.open(cns_out_fpath, "rb")
    else:
        file_handle = open(cns_out_fpath, "rb")
    # Read the file
    try:
        _find_cns_errors(file_handle, KNOWN_ERRORS, filepath=cns_out_fpath)
    except KnownCNSError as err:
        return err
    else:
        # return the cause
        return KnownCNSError(
            "An unfortunate CNS error occured at exection time...",
            f"Manually check the file `{cns_out_fpath}` to understand why!",
            cns_out_fpath,
            )


def _find_cns_errors(
        file_handle: Union[gzip.GzipFile, BufferedReader],
        known_errors: dict[str, str],
        chunk_size: int = 4096,
        filepath: FilePath = "",
        ) -> None:
    """Backward reading and detect first known CNS error in file.

    Parameters
    ----------
    file_handle: Union[gzip.GzipFile, BufferedReader]
        An opened file in read bytes mode.
    known_errors : dict[str, str]
        Dict of known errors and their hints
    chunk_size : int, optional
        Check size (in bytes) to read the file backwards, by default 4096
    filepath : FilePath -> Union[str, Path]
        Path to the cns.cnserr file currently checked.

    Raises
    ------
    KnownCNSError
        An exception for known CNS errors, with its hint on how to solve it!
    """
    # Find file size
    file_handle.seek(0, 2)
    size = file_handle.tell()
    buffer = b''
    # Set the number of lines to parse
    # initiated with high value to read all lines the first time
    # updated to the number of lines that were already parsed
    # after each chunk iteration, so we do not parse the same lines
    parsed_lines = 99999
    for i in range(size - 1, -1, -chunk_size):
        # Go to location in file
        file_handle.seek(max(i - chunk_size, 0))
        # Read next chunk
        chunk = file_handle.read(min(chunk_size, i + 1))
        # Increment buffer
        buffer = chunk + buffer
        lines = buffer.split(b'\n')
        # Read lines
        for line in reversed(lines[-len(lines):parsed_lines]):
            decoded_line = line.decode('utf-8', errors='replace')
            # Loop over known errors
            for error_string, hint in known_errors.items():
                # Check if this error is known
                if error_string in decoded_line:
                    # return the cause
                    raise KnownCNSError(
                        error_string,
                        hint,
                        filepath,
                        )
        # Update number of parsed lines so we do not check them again
        parsed_lines = -len(lines)


def find_all_cns_errors(
        directory_path: FilePath,
        ) -> dict[str, dict[str, Union[int, KnownCNSError]]]:
    """Find all errors in a directory.

    Parameters
    ----------
    directory_path : FilePath
        Path to the directory to be checked.

    Returns
    -------
    all_errors : dict[str, dict[str, Union[list[FilePath], KnownCNSError]]]
        Dictionary containing all errors found in this directory.
    """
    all_errors: dict[str, dict[str, Union[int, KnownCNSError]]] = {}
    # Gather list of all `.cnserr` and `.cnserr.gz` files present in directory
    all_cns_out_files = list(Path(directory_path).glob("*.cnserr.gz"))
    all_cns_out_files += list(Path(directory_path).glob("*.cnserr"))
    # Loop over all .cnserr files
    for fpath in all_cns_out_files:
        # Try to dectect an error
        if (detected_error := find_cns_errors(fpath)):
            # Hold data if an error is present in that file
            error_type = all_errors.setdefault(
                detected_error.cns_error,
                {"files": [], "error": detected_error}
                )
            error_type["files"].append(detected_error.filepath)
    return all_errors
