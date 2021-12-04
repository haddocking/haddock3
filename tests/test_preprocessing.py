"""Test preprocessing operations."""
import logging
import os
from pathlib import Path

import pytest

from haddock.gear import preprocessing as pp

from . import broken_pdb, data_folder, good_pdb


def test_all(caplog):
    """."""
    caplog.set_level(logging.WARNING)

    pp.process_pdb_files(broken_pdb)

    output = Path(
        broken_pdb.parent,
        broken_pdb.stem + '_processed').with_suffix(broken_pdb.suffix)

    lines_results = output.read_text().split('\n')
    lines_expected = good_pdb.read_text().split('\n')

    # asserts the number of WARNING messages is appropriate
    # see: https://docs.pytest.org/en/6.2.x/logging.html
    assert len(caplog.text.split(os.linesep)) > 0

    assert lines_results[-1] == lines_expected[-1]
    assert len(lines_results) == len(lines_expected)
    for i in range(len(lines_results)):
        assert lines_results[i] == lines_expected[i]
