"""Test client copy."""
import shutil
from pathlib import Path

from haddock.clis.cli_cp import main

from . import tests_path


def test_main():
    """Test main function of haddock3-copy CLI."""
    run1 = Path(tests_path, "clis", "hd3_copy", "run1")
    run2 = Path(tests_path, "clis", "hd3_copy", "run2")
    run2_ref = Path(tests_path, "clis", "hd3_copy", "run2_ref")
    main(run1, [0, 2, 3], run2)
    compare_files(run2, run2_ref)
    # easy way to make sure no file is left behind
    compare_files(run2_ref, run2)
    shutil.rmtree(run2)


def compare_files(folder1, folder2):
    for file_folder in folder1.iterdir():
        file_folder_2 = Path(folder2, file_folder.name)
        if file_folder.is_file():
            text1 = file_folder.read_text()
            text2 = file_folder_2.read_text()
            assert text1 == text2
        if file_folder.is_dir():
            compare_files(file_folder, file_folder_2)
