"""Test preprocessing client."""

import os
from pathlib import Path

import pytest
import tempfile
import shutil
from haddock.clis.cli_pp import ap, main

from . import data_folder


@pytest.fixture(name="broken_pdb")
def fixture_broken_pdb():
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)

        src = Path(data_folder, "broken.pdb")
        dst = Path(tempdir, "broken.pdb")
        shutil.copy(src, dst)

        yield dst


@pytest.fixture(name="corrected_pdb")
def fixture_corrected_pdb():
    with (
        tempfile.NamedTemporaryFile(delete=False) as temp_file,
        open(Path(data_folder, "corrected.pdb"), "rb") as source_file,
    ):
        shutil.copyfileobj(source_file, temp_file)
        yield Path(temp_file.name)


def test_pp_cli(broken_pdb, corrected_pdb):
    """
    Test `haddock-pp` client with broken PDB.

    This is the same PDB used to test the preprocessing library.
    """
    main(broken_pdb)
    output = Path(f"{Path(broken_pdb).stem}_processed.pdb")
    assert output.exists()

    result = output.read_text().strip(os.linesep)
    expected = corrected_pdb.read_text().strip(os.linesep)
    assert result == expected


@pytest.mark.parametrize(
    "arg,key,value",
    [
        ("-d", "dry", True),
        ("--dry", "dry", True),
        ("-t file1.top file2.top", "topfile", ["file1.top", "file2.top"]),
        (
            "--topfile file1.top file2.top",
            "topfile",
            ["file1.top", "file2.top"],
        ),  # noqa: E501
        ("-s somesuffix", "suffix", "somesuffix"),
        ("--suffix somesuffix", "suffix", "somesuffix"),
        ("-odir somedir", "output_directory", Path("somedir")),
        ("--output-directory somedir", "output_directory", Path("somedir")),
    ],
)
def test_cli_args(arg, key, value):
    """Test adding arguments to client."""
    cmd = ap.parse_args(("some-file.pdb " + arg).split())
    assert value == vars(cmd)[key]
