"""Unitests related to haddock.gear.known_cns_errors.py."""

import pytest
import tempfile

from pathlib import Path

from haddock.gear.known_cns_errors import (
    KNOWN_ERRORS,
    find_cns_errors,
    find_all_cns_errors,
    )
from haddock.core.exceptions import KnownCNSError


def test_find_cns_errors():
    with tempfile.TemporaryDirectory("errored_module") as tmp:
        errored_filepath = Path(tmp, "with_error_cns.out")
        for error, message in KNOWN_ERRORS.items():
            # Write error in a file
            errored_filepath.write_text(error)
            # Detect it
            detected_error = find_cns_errors(errored_filepath)
            # Check error was detected
            assert type(detected_error) == KnownCNSError
            # Check hint message contained in error string representation
            assert message in str(detected_error)


def test_find_all_cns_errors():
    with tempfile.TemporaryDirectory("errored_module") as tmp:
        for i, error in enumerate(KNOWN_ERRORS.keys()):
            # Create two files with same error
            for j in range(2):
                errored_filepath = Path(tmp, f"with_error_cns_{i}_{j}.out")
                # Write error in a file
                errored_filepath.write_text(error)
        # Detect all errors within a directory
        all_errors = find_all_cns_errors(tmp)
        # Check not empty
        assert all_errors != {}
        # Check that all keys are present
        assert len(all_errors.keys()) == len(KNOWN_ERRORS.keys())
        # Loop over detected errors
        for cns_error in all_errors.keys():
            error = all_errors[cns_error]
            # Check that we detected both files
            assert error["count"] == 2
            # Check that error hint is well reported
            assert KNOWN_ERRORS[cns_error] in str(error["error"])
