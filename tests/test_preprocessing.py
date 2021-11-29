"""Test preprocessing operations."""
from pathlib import Path

import pytest

from haddock.gear import preprocessing as pp

from . import broken_pdb, data_folder, good_pdb


def test_all():
    """."""
    pp.process_pdb_files(broken_pdb)

    output = Path(
        broken_pdb.parent,
        broken_pdb.stem + '_processed').with_suffix(broken_pdb.suffix)

    lines_results = output.read_text().split('\n')
    lines_expected = good_pdb.read_text().split('\n')

    assert len(lines_results) == len(lines_expected)
    for i in range(len(lines_results)):
        assert lines_results[i] == lines_expected[i]
