"""Integration tests related to haddock.gear.known_cns_errors.py."""

import gzip
import pytest
import tempfile
import random

from os import linesep
from pathlib import Path
from string import ascii_letters

from haddock.gear.known_cns_errors import KNOWN_ERRORS
from haddock.libs.libontology import PDBFile
from haddock.modules.sampling.rigidbody import (
    DEFAULT_CONFIG as DEFAULT_RIGIDBODY_CONFIG,
    HaddockModule as RigidbodyModule
    )


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
    """Generate directory full of CNS.cnserr file with errors."""
    with tempfile.TemporaryDirectory("moduleoutputs") as tmp:
        for i, error in enumerate(KNOWN_ERRORS.keys()):
            # Generate an error string in the middle of the file
            error_text = gen_random_text + error + gen_random_text
            # Create two files with same error
            for j in range(1, 3):
                errored_filepath = Path(tmp, f"with_error_cns_{i}_{j}.cnserr")
                # Write error in a file
                errored_filepath.write_text(error_text)
            # Create two compressed files with same error
            for j in range(1, 3):
                errored_gz_file = Path(tmp, f"with_error_cns_{i}_{j}.cnserr.gz")
                # Write error in a file
                with gzip.open(errored_gz_file, mode="wb") as gout:
                    gout.write(bytes(error_text, encoding="utf-8"))
        yield tmp


@pytest.fixture
def rigidbody_module_with_cns_errors(gen_fake_cns_errors):
    """Generate a failed rigidbody module with CNS errors."""
    rigidbody = RigidbodyModule(
        order=1,
        path=Path(gen_fake_cns_errors),
        initial_params=DEFAULT_RIGIDBODY_CONFIG,
        )
    # Generate 9 filepath that were not created
    rigidbody.output_models = [
        PDBFile(Path(gen_fake_cns_errors, f"not_generated_output_{i}.pdb"))
        for i in range(1, 10)
        ]
    yield rigidbody


@pytest.fixture
def rigidbody_module_without_cns_errors():
    """Generate a failed rigidbody module without CNS errors."""
    with tempfile.TemporaryDirectory("moduleoutputs") as tmp:
        rigidbody = RigidbodyModule(
            order=1,
            path=Path(tmp),
            initial_params=DEFAULT_RIGIDBODY_CONFIG,
            )
        # Generate 9 filepath that were not created
        rigidbody.output_models = [
            PDBFile(Path(tmp, f"not_generated_output_{i}.pdb"))
            for i in range(1, 10)
            ]
        yield rigidbody


class MockPreviousIO:
    """Mock proviousIO function."""

    def __init__(self, path):
        self.path = path
        self.output = []


def test_detection_when_faulty(rigidbody_module_with_cns_errors):
    """Test failure of run and detection of CNS errors."""
    rigidbody_module_with_cns_errors.previous_io = MockPreviousIO(
        rigidbody_module_with_cns_errors.path
        )
    # Check that the run will fail
    with pytest.raises(RuntimeError) as error_info:
        rigidbody_module_with_cns_errors.export_io_models()
    # Get final error string
    string_error = str(error_info.value)
    # Loop over known errors
    for cns_error_string, user_hint in KNOWN_ERRORS.items():
        # Check it was detected
        assert cns_error_string in string_error
        # Check user hint is present in error message
        assert user_hint in string_error


def test_undetected_when_faulty(rigidbody_module_without_cns_errors):
    """Test failure of run and undetection of CNS errors."""
    rigidbody_module_without_cns_errors.previous_io = MockPreviousIO(
        rigidbody_module_without_cns_errors.path
        )
    # Check that the run will fail
    with pytest.raises(RuntimeError) as error_info:
        rigidbody_module_without_cns_errors.export_io_models()
    # Get final error string
    string_error = str(error_info.value)
    # Loop over known errors
    for cns_error_string, user_hint in KNOWN_ERRORS.items():
        # Check it was NOT detected
        assert cns_error_string not in string_error
        # Check user hint NOT is present in error message
        assert user_hint not in string_error
