"""Test preprocessing client."""
import os
from pathlib import Path

import pytest

from haddock.clis.cli_pp import ap, main

from . import broken_pdb, corrected_pdb


def test_pp_cli():
    """
    Test `haddock-pp` client with broken PDB.

    This is the same PDB used to test the preprocessing library.
    """
    main(broken_pdb)
    output = Path(Path.cwd(), 'broken_processed.pdb')
    assert output.exists()

    result = output.read_text().strip(os.linesep)
    expected = corrected_pdb.read_text().strip(os.linesep)
    assert result == expected

    output.unlink()


@pytest.mark.parametrize(
    "arg,key,value",
    [
        ('-d', "dry", True),
        ('--dry', "dry", True),
        ('-t file1.top file2.top', "topfile", ['file1.top', 'file2.top']),
        ('--topfile file1.top file2.top', "topfile", ['file1.top', 'file2.top']),  # noqa: E501
        ('-s somesuffix', 'suffix', 'somesuffix'),
        ('--suffix somesuffix', 'suffix', 'somesuffix'),
        ('-odir somedir', 'output_directory', Path('somedir')),
        ('--output-directory somedir', 'output_directory', Path('somedir')),
        ]
    )
def test_cli_args(arg, key, value):
    """Test adding arguments to client."""
    cmd = ap.parse_args(('some-file.pdb ' + arg).split())
    assert value == vars(cmd)[key]
