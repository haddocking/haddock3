"""Unitests related to haddock.gear.known_cns_errors.py."""

import gzip
import pytest
import tempfile
import random

from os import linesep
from pathlib import Path
from string import ascii_letters

from haddock.gear.known_cns_errors import (
    KNOWN_ERRORS,
    find_cns_errors,
    find_all_cns_errors,
    )
from haddock.core.exceptions import KnownCNSError


@pytest.fixture
def gen_random_text():
    """Generate some random text."""
    textline = "".join([random.choice(ascii_letters) for _ in range(80)])
    text = ""
    for _ in range(500):
        text += f"{textline}{linesep}"
    yield text


@pytest.fixture
def gen_fake_cns_errors(gen_random_text):
    """Generate directory full of CNS.out file with errors."""
    with tempfile.TemporaryDirectory("moduleoutputs") as tmp:
        for i, error in enumerate(KNOWN_ERRORS.keys()):
            error_text = gen_random_text + error + gen_random_text
            # Create two files with same error
            for j in range(1, 3):
                errored_filepath = Path(tmp, f"errored_cns_{i}_{j}.cnserr")
                # Write error in a file
                errored_filepath.write_text(error_text)
            # Create two compressed files with same error
            for j in range(1, 3):
                errored_gz_file = Path(tmp, f"errored_cns_{i}_{j}.cnserr.gz")
                # Write error in a file
                with gzip.open(errored_gz_file, mode="wb") as gout:
                    gout.write(bytes(error_text, encoding="utf-8"))
        yield tmp


def test_find_cns_errors(gen_random_text):
    """Test detection of error in a file."""
    with tempfile.TemporaryDirectory("errored_module") as tmp:
        errored_filepath = Path(tmp, "errored_cns.cnserr")
        for error, message in KNOWN_ERRORS.items():
            # Write error in a file
            error_text = gen_random_text + error + gen_random_text
            errored_filepath.write_text(error_text)
            # Detect it
            detected_error = find_cns_errors(errored_filepath)
            # Check error was detected
            assert type(detected_error) == KnownCNSError
            # Check hint message contained in error string representation
            assert message in str(detected_error)


def test_find_all_cns_errors(gen_fake_cns_errors):
    """Test detection of multiple errors in a directory."""
    # Detect all errors within a directory
    all_errors = find_all_cns_errors(gen_fake_cns_errors)
    # Check not empty
    assert all_errors != {}
    # Check that all keys are present
    assert len(all_errors.keys()) == len(KNOWN_ERRORS.keys())
    # Loop over detected errors
    for cns_error in all_errors.keys():
        error = all_errors[cns_error]
        # Check that we detected both files (for both .cnserr and .cnserr.gz)
        assert len(error["files"]) == 4  # 2 * 2
        # Check that error hint is well reported
        assert KNOWN_ERRORS[cns_error] in str(error["error"])
