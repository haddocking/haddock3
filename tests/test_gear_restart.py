"""Test gear.restart_run."""
import argparse
import os
import pytest
import tempfile

from haddock.gear import (
    restart_run,
    preprocess_restart_from,
    ANA_DIR,
    TRACEBACK_FOLDER,
)


def test_has_help():
    """Assert module has _help_cli variable."""
    assert restart_run._help_cli


@pytest.mark.parametrize(
    'i,expected',
    [
        ('0', 0),
        ('1', 1),
        ('57', 57),
        (100, 100),
        ]
    )
def test_non_neg_arg(i, expected):
    """Test non negative arg type."""
    r = restart_run._arg_non_neg_int(i)
    assert r == expected


@pytest.mark.parametrize(
    'i,expected',
    [
        ('0', 0),
        ('1', 1),
        ('57', 57),
        (100, 100),
        ]
    )
def test_restart_cli(i, expected):
    """Test non negative arg type."""
    ap = argparse.ArgumentParser()
    restart_run.add_restart_arg(ap)
    cmd = ap.parse_args(f'--restart {i}'.split())
    assert cmd.restart == expected


@pytest.mark.parametrize(
    'n',
    (-1, -10, '-1230', -50000),
    )
def test_arg_non_neg_error(n):
    """Test non-negative int."""
    with pytest.raises(argparse.ArgumentTypeError):
        restart_run._arg_non_neg_int(n)


@pytest.mark.parametrize(
    'n',
    (-1, -10, '-1230', -50000),
    )
def test_restart_cli_error(n):
    """Test --restart number with error."""
    ap = argparse.ArgumentParser()
    restart_run.add_restart_arg(ap)
    with pytest.raises(SystemExit) as exit:
        ap.parse_args(f'--restart {n}'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 2


def test_preprocess_restart_from():
    """Test removal of downstream directories."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Build mimic of previous run directories
        module_dirs = ["topoaa", "rigidbody", "caprieval"]
        for module_index, module_name in module_dirs:
            os.makedirs(f"{tmpdir}/data/{module_index}_{module_name}")
            os.makedirs(f"{tmpdir}/{module_index}_{module_name}")
            if module_name == "caprieval":
                os.makedirs(
                    f"{tmpdir}/{ANA_DIR}/{module_index}_{module_name}_analysis"
                    )
        os.makedirs(f"{tmpdir}/{TRACEBACK_FOLDER}")
        
        # Test restart preprocessing function
        preprocess_restart_from(tmpdir, 2)
        
        # Verify they were removed
        expected_to_be_removed = [
            f"{tmpdir}/{TRACEBACK_FOLDER}",
            f"{tmpdir}/{ANA_DIR}/2_caprieval_analysis",
            f"{tmpdir}/2_caprieval",
        ]
        expected_to_stay_there = [
            f"{tmpdir}/1_rigidbody",
            f"{tmpdir}/data/1_rigidbody",
            f"{tmpdir}/0_topoaa",
            f"{tmpdir}/data/0_topoaa",
        ]
        for should_not_exist in expected_to_be_removed:
            assert not os.path.exists(should_not_exist)
        for should_exist in expected_to_stay_there:
            assert os.path.exists(should_exist)
