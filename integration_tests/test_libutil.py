from pathlib import Path

import pytest

from haddock.libs.libutil import get_cns_executable

from . import is_linux_x86_64, is_not_linux_x86_64


def test_finds_default_cns_executable():
    """Test that the function finds CNS in the default location."""
    cns_exec, cns_exec_linux = get_cns_executable()

    assert cns_exec.exists(), f"CNS executable not found at {cns_exec}"
    assert isinstance(cns_exec, Path)
    assert isinstance(cns_exec_linux, Path)


@is_linux_x86_64
def test_cns_executable_is_the_same():
    """Test that the function finds CNS in the default location."""
    cns_exec, cns_exec_linux = get_cns_executable()

    assert cns_exec == cns_exec_linux


def test_cns_executable_is_file():
    """Test that the CNS executable is actually a file."""
    cns_exec, _ = get_cns_executable()

    assert cns_exec.is_file(), f"{cns_exec} is not a file"


@is_linux_x86_64
def test_linux_x86_64_returns_same_executable():
    """On x86_64 Linux, both executables should be the same."""
    cns_exec, cns_exec_linux = get_cns_executable()

    assert cns_exec == cns_exec_linux


@is_not_linux_x86_64
def test_non_linux_looks_for_linux_variant():
    """On non-Linux systems, should look for _linux variant."""
    cns_exec, cns_exec_linux = get_cns_executable()

    assert cns_exec != cns_exec_linux
    assert cns_exec_linux.name == f"{cns_exec.name}_linux"


def test_returns_tuple_of_paths():
    """Test that function returns a tuple of two Path objects."""
    result = get_cns_executable()

    assert isinstance(result, tuple)
    assert len(result) == 2
    assert all(isinstance(p, Path) for p in result)
