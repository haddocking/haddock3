"""This file contains functions to detect known/common CNS errors.

Inspired from:
https://github.com/haddocking/haddock25/blob/main/tools/check-error-messages.sh
"""

from haddock.core.exceptions import KnownCNSError
from haddock.core.typing import Optional

# Dictionary of known errors
# as key:    How to catch it in the cns.out
# as value:  Message to user
KNOWN_ERRORS = {
    "CHAIN LENGHT FOR SYMMETRY RESTRAINTS DO NOT MATCH": (
        "Missmatch between chain length for symmetry restraints. "
        "Check your input molecules and symmetry restraints."
    ),
    "NCS-restraints error encountered: Improperly defined non-crystallographic symmetry": (
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
        "Make sure to filter those for solvent accessibility."
        ),
    "SELRPN error encountered: parsing error": (
        "Check your restraint files."
    ),
    "PARSER error encountered: Encountered too many parsing errors": (
        "Encountered too many parsing errors."
        ),
    "XMREAD error encountered:  sectioning of map incompatible with resolution": (
        "Check your EM map resolution and sectioning"
        ),
    "ALLHP error encountered: not enough memory available": (
        "Too many distance restraints defined. "
        "Try to reduce this number by checking your definition of active and "
        "passive residues. "
        "Make sure to filter those for solvent accessibility. "
        "Try to decrease the size of your system where possible."
        ),
    "error encountered: missing SCATter definition for SELEcted atoms": (
        "Unsupported atoms/molecules for cryo-EM restraints"
        ),
    "ROTMAT error encountered: rotation vector has zero length": (
        "Check your input parameters and restraints"
        "Possibly try turning off the sampling of 180 degrees rotattion"
        )
    }  # noqa : E501


def find_cns_errors(cns_out_fpath: str) -> Optional[KnownCNSError]:
    """Detect if a known CNS error is in a cns.out file.

    Parameters
    ----------
    cns_out_fpath : str
        Path to the cns.out file to check.

    Returns
    -------
    Optional[KnownCNSError]
        An exception for known CNS errors, with its hint on how to solve it!
    """
    try:
        _find_cns_errors(cns_out_fpath, KNOWN_ERRORS)
    except KnownCNSError as err:
        return err
    else:
        return None


def _find_cns_errors(
        cns_out_fpath: str,
        known_errors: dict[str, str],
        chunk_size: int = 4096,
        ) -> None:
    """Backward reading and detect first known CNS error in file.

    Parameters
    ----------
    cns_out_fpath : str
        Path to the cns.out file to check.
    known_errors : dict[str, str]
        Dict of known errors and their hints
    chunk_size : int, optional
        Check size (in bytes) to read the file backwards, by default 4096

    Raises
    ------
    KnownCNSError
        An exception for known CNS errors, with its hint on how to solve it!
    """
    # Read file
    with open(cns_out_fpath, 'rb') as file:
        # Find file size
        file.seek(0, 2)
        size = file.tell()
        buffer = b''
        parsed_lines = 9999
        for i in range(size - 1, -1, -chunk_size):
            # Go to location in file
            file.seek(max(i - chunk_size, 0))
            # Read next chunk
            chunk = file.read(min(chunk_size, i + 1))
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
                            )
            # Update number of parsed lines so we do not check them again
            parsed_lines = -len(lines)
