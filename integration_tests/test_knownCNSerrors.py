"""Integration tests related to haddock.gear.known_cns_errors.py."""

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
    """Generate directory full of CNS.out file with errors."""
    with tempfile.TemporaryDirectory("moduleoutputs") as tmp:
        for i, error in enumerate(KNOWN_ERRORS.keys()):
            error_text = gen_random_text + error + gen_random_text
            # Create two files with same error
            for j in range(2):
                errored_filepath = Path(tmp, f"with_error_cns_{i}_{j}.out")
                # Write error in a file
                errored_filepath.write_text(error_text)
        yield tmp


@pytest.fixture
def rigidbody_module(gen_fake_cns_errors):
    rigidbody = RigidbodyModule(
        order=1,
        path=Path(gen_fake_cns_errors),
        initial_params=DEFAULT_RIGIDBODY_CONFIG,
        )
    rigidbody.output_models = [
        PDBFile(Path(gen_fake_cns_errors, f"none_generated_output_{i}.pdb"))
        for i in range(10)
        ]
    yield rigidbody


class MockPreviousIO:
    def __init__(self, path):
        self.path = path
        self.output = []


def test_detection_when_faulty(rigidbody_module):
    rigidbody_module.previous_io = MockPreviousIO(rigidbody_module.path)
    # Check that the run will fail
    with pytest.raises(RuntimeError) as error_info:
        rigidbody_module.export_io_models()
    # Get final error string
    string_error = str(error_info.value)
    # Loop over known errors
    for cns_error_string, user_hint in KNOWN_ERRORS.items():
        # Check it was detected
        assert cns_error_string in string_error
        # Check user hint is present in error message
        assert user_hint in string_error