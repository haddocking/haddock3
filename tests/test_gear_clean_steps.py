"""Test clean steps."""
import shutil
from pathlib import Path

from haddock.gear.clean_steps import (
    clean_output,
    unpack_compressed_and_archived_files,
    update_unpacked_names,
    )

from . import clean_steps_folder


def test_clean_output():
    """Test correct clean and unpack functions."""
    # defines the dir to compress
    outdir = Path(clean_steps_folder, 'run1c')
    # remove it just in case it remained in the testing environment
    shutil.rmtree(outdir, ignore_errors=True)

    # copies the directory to avoid messing with git stack
    shutil.copytree(
        Path(clean_steps_folder, 'run1'),
        outdir,
        )

    # packs each folder separately
    clean_output(Path(outdir, "0_topoaa"))
    clean_output(Path(outdir, "1_rigidbody"))
    clean_output(Path(outdir, "2_clustfcc"))

    # these are the files that SHOULD exist after packing
    files_that_should_exist = [
        Path("0_topoaa", "params.cfg"),
        Path("0_topoaa", "structure_1.inp.gz"),
        Path("0_topoaa", "structure_1.out.gz"),
        Path("0_topoaa", "structure_1.pdb.gz"),
        Path("0_topoaa", "structure_1.psf.gz"),

        Path("1_rigidbody", "params.cfg"),
        Path("1_rigidbody", "io.json"),
        Path("1_rigidbody", "structure_1.inp.gz"),
        Path("1_rigidbody", "structure_1.out.gz"),
        Path("1_rigidbody", "seed.tgz"),
        Path("1_rigidbody", "structure_1.pdb.gz"),
        Path("1_rigidbody", "structure_2.pdb.gz"),

        Path("2_clustfcc", "con.tgz"),
        ]

    # assert if the files actually exist in the `outdir` folder
    for file_ in files_that_should_exist:
        assert Path(outdir, file_).exists()

    # these are the files that should NOT exist after packing
    files_that_should_not_exist = [
        Path("0_topoaa", "structure_1.inp"),
        Path("0_topoaa", "structure_1.out"),
        Path("0_topoaa", "structure_1.pdb"),
        Path("0_topoaa", "structure_1.psf"),

        Path("1_rigidbody", "structure_1.inp"),
        Path("1_rigidbody", "structure_1.out"),
        Path("1_rigidbody", "structure_1.pdb"),
        Path("1_rigidbody", "structure_1.seed"),
        Path("1_rigidbody", "structure_2.inp"),
        Path("1_rigidbody", "structure_2.out"),
        Path("1_rigidbody", "structure_2.pdb"),
        Path("1_rigidbody", "structure_2.seed"),

        Path("2_clustfcc", "structure_2.con"),
        Path("2_clustfcc", "structure_2.con"),
        ]

    # asserts the files do NOT exist in the `outdir` folder
    for file_ in files_that_should_not_exist:
        assert not Path(outdir, file_).exists()

    # now, unpacks the `outdir`. In other words, reverts the previous
    # clean_output operation
    unpack_compressed_and_archived_files([
        Path(outdir, "0_topoaa"),
        Path(outdir, "1_rigidbody"),
        Path(outdir, "2_clustfcc"),
        ])

    # these are the files that SHOULD exist after unpacking
    files_that_should_exist = [
        Path("0_topoaa", "params.cfg"),
        # Path("0_topoaa", "structure_1.inp"),
        # Path("0_topoaa", "structure_1.out"),
        Path("0_topoaa", "structure_1.pdb"),
        Path("0_topoaa", "structure_1.psf"),
        Path("0_topoaa", "structure_1.inp.gz"),
        Path("0_topoaa", "structure_1.out.gz"),

        Path("1_rigidbody", "params.cfg"),
        Path("1_rigidbody", "io.json"),
        # Path("1_rigidbody", "structure_1.inp"),
        # Path("1_rigidbody", "structure_1.out"),
        Path("1_rigidbody", "structure_1.pdb"),
        Path("1_rigidbody", "structure_1.seed"),
        Path("1_rigidbody", "structure_2.pdb"),
        Path("1_rigidbody", "structure_2.seed"),

        Path("2_clustfcc", "structure_1.con"),
        Path("2_clustfcc", "structure_2.con"),
        ]

    # these are the files that should NOT exist after unpacking
    files_that_should_not_exist = [
        Path("0_topoaa", "structure_1.inp"),
        Path("0_topoaa", "structure_1.out"),
        Path("0_topoaa", "structure_1.pdb.gz"),
        Path("0_topoaa", "structure_1.psf.gz"),

        Path("1_rigidbody", "structure_1.inp"),
        Path("1_rigidbody", "structure_1.out"),
        Path("1_rigidbody", "seed.tgz"),
        Path("1_rigidbody", "structure_1.pdb.gz"),
        Path("1_rigidbody", "structure_2.pdb.gz"),

        Path("2_clustfcc", "con.tgz"),
        ]

    # assert files actually exist
    for file_ in files_that_should_exist:
        assert Path(outdir, file_).exists()

    # assert files do not exist
    for file_ in files_that_should_not_exist:
        assert not Path(outdir, file_).exists()

    # removes the dir as it is no longer needed
    shutil.rmtree(outdir)


def test_clean_output_dec_all():
    """Test correct clean and unpack functions."""
    # defines the dir to compress
    outdir = Path(clean_steps_folder, 'run1c')
    # remove it just in case it remained in the testing environment
    shutil.rmtree(outdir, ignore_errors=True)

    # copies the directory to avoid messing with git stack
    shutil.copytree(
        Path(clean_steps_folder, 'run1'),
        outdir,
        )

    # packs each folder separately
    clean_output(Path(outdir, "0_topoaa"))
    clean_output(Path(outdir, "1_rigidbody"))
    clean_output(Path(outdir, "2_clustfcc"))

    # these are the files that SHOULD exist after packing
    files_that_should_exist = [
        Path("0_topoaa", "params.cfg"),
        Path("0_topoaa", "structure_1.inp.gz"),
        Path("0_topoaa", "structure_1.out.gz"),
        Path("0_topoaa", "structure_1.pdb.gz"),
        Path("0_topoaa", "structure_1.psf.gz"),

        Path("1_rigidbody", "params.cfg"),
        Path("1_rigidbody", "io.json"),
        Path("1_rigidbody", "structure_1.inp.gz"),
        Path("1_rigidbody", "structure_1.out.gz"),
        Path("1_rigidbody", "seed.tgz"),
        Path("1_rigidbody", "structure_1.pdb.gz"),
        Path("1_rigidbody", "structure_2.pdb.gz"),

        Path("2_clustfcc", "con.tgz"),
        ]

    # assert if the files actually exist in the `outdir` folder
    for file_ in files_that_should_exist:
        assert Path(outdir, file_).exists()

    # these are the files that should NOT exist after packing
    files_that_should_not_exist = [
        Path("0_topoaa", "structure_1.inp"),
        Path("0_topoaa", "structure_1.out"),
        Path("0_topoaa", "structure_1.pdb"),
        Path("0_topoaa", "structure_1.psf"),

        Path("1_rigidbody", "structure_1.inp"),
        Path("1_rigidbody", "structure_1.out"),
        Path("1_rigidbody", "structure_1.pdb"),
        Path("1_rigidbody", "structure_1.seed"),
        Path("1_rigidbody", "structure_2.inp"),
        Path("1_rigidbody", "structure_2.out"),
        Path("1_rigidbody", "structure_2.pdb"),
        Path("1_rigidbody", "structure_2.seed"),

        Path("2_clustfcc", "structure_2.con"),
        Path("2_clustfcc", "structure_2.con"),
        ]

    # asserts the files do NOT exist in the `outdir` folder
    for file_ in files_that_should_not_exist:
        assert not Path(outdir, file_).exists()

    # now, unpacks the `outdir`. In other words, reverts the previous
    # clean_output operation
    unpack_compressed_and_archived_files([
        Path(outdir, "0_topoaa"),
        Path(outdir, "1_rigidbody"),
        Path(outdir, "2_clustfcc"),
        ],
        dec_all=True)

    # these are the files that SHOULD exist after FULL unpacking
    files_that_should_exist = [
        Path("0_topoaa", "params.cfg"),
        Path("0_topoaa", "structure_1.inp"),
        Path("0_topoaa", "structure_1.out"),
        Path("0_topoaa", "structure_1.pdb"),
        Path("0_topoaa", "structure_1.psf"),

        Path("1_rigidbody", "params.cfg"),
        Path("1_rigidbody", "io.json"),
        Path("1_rigidbody", "structure_1.inp"),
        Path("1_rigidbody", "structure_1.out"),
        Path("1_rigidbody", "structure_1.pdb"),
        Path("1_rigidbody", "structure_1.seed"),
        Path("1_rigidbody", "structure_2.pdb"),
        Path("1_rigidbody", "structure_2.seed"),

        Path("2_clustfcc", "structure_1.con"),
        Path("2_clustfcc", "structure_2.con"),
        ]

    # these are the files that should NOT exist after unpacking
    files_that_should_not_exist = [
        Path("0_topoaa", "structure_1.inp.gz"),
        Path("0_topoaa", "structure_1.out.gz"),
        Path("0_topoaa", "structure_1.pdb.gz"),
        Path("0_topoaa", "structure_1.psf.gz"),

        Path("1_rigidbody", "structure_1.inp.gz"),
        Path("1_rigidbody", "structure_1.out.gz"),
        Path("1_rigidbody", "seed.tgz"),
        Path("1_rigidbody", "structure_1.pdb.gz"),
        Path("1_rigidbody", "structure_2.pdb.gz"),

        Path("2_clustfcc", "con.tgz"),
        ]

    # assert files actually exist
    for file_ in files_that_should_exist:
        assert Path(outdir, file_).exists()

    # assert files do not exist
    for file_ in files_that_should_not_exist:
        assert not Path(outdir, file_).exists()

    # removes the dir as it is no longer needed
    shutil.rmtree(outdir)


def test_update_unpacked_names():
    original = ['0_topoaa', Path('4_flexref'), '5_seletopclusts']
    prev = ['0_topoaa', 'run_dir/4_flexref', '5_seletopclusts']
    new = ['run_dir/0_topoaa', '1_flexref', 'run_dir/2_seletopclusts']
    update_unpacked_names(prev, new, original)
    assert original == ['0_topoaa', Path('1_flexref'), '2_seletopclusts']
