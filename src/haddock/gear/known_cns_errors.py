"""This file contains functions to detect known/common CNS errors.

Inspired from:
https://github.com/haddocking/haddock25/blob/main/tools/check-error-messages.sh
"""

from haddock.core.exceptions import KnownCNSError
from haddock.core.typing import Optional

# Dictionary of known errors
# as key:    How to catch it in the cns.out
# as value:  Message to user
known_errors = {
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
    if (detected_error := _find_cns_errors(cns_out_fpath, known_errors)):
        return KnownCNSError(known_errors[detected_error])


def _find_cns_errors(
        cns_out_fpath: str,
        known_errors: dict[str, str],
        chunk_size: int = 4096,
        ) -> Optional[str]:
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
                for error_string in known_errors.keys():
                    # Check if this error is known
                    if error_string in decoded_line:
                        # return the cause
                        return error_string
            # Update number of parsed lines so we do not check them again
            parsed_lines = -len(lines)
