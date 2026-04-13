import os
from importlib.resources import files
from pathlib import Path

import pytest

import haddock
from haddock.libs.libutil import get_cns_executable, get_prodrg_exec

from . import is_linux_x86_64, is_not_linux_x86_64


def test_finds_default_cns_executable():
    """Test that the function finds CNS in the default location."""
    cns_exec, cns_exec_linux = get_cns_executable()

    assert cns_exec.exists(), f"CNS executable not found at {cns_exec}"
    assert isinstance(cns_exec, Path)
    assert isinstance(cns_exec_linux, Path)


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
    assert cns_exec_linux.name == "x86_64-linux.bin"


def test_cns_binary_is_executable():
    """Make sure the binaries are executable."""
    cns_exec_dir = Path(files(haddock).joinpath("cns/bin"))  # type: ignore
    for f in cns_exec_dir.glob("*"):
        assert os.access(f, os.X_OK)


def test_returns_tuple_of_paths():
    """Test that function returns a tuple of two Path objects."""
    result = get_cns_executable()

    assert isinstance(result, tuple)
    assert len(result) == 2
    assert all(isinstance(p, Path) for p in result)


def test_get_prodrg_exec_returns_tuple():
    """Test that get_prodrg_exec returns a tuple of length 2."""
    result = get_prodrg_exec()

    assert isinstance(result, tuple)
    assert len(result) == 2


def test_prodrg_exec_is_file_when_available():
    """If prodrg is available, the executable should be a real file."""
    prodrg_exec, _ = get_prodrg_exec()

    if prodrg_exec is not None:
        assert prodrg_exec.is_file(), f"{prodrg_exec} is not a file"


def test_prodrg_param_is_file_when_available():
    """If prodrg is available, the param file should be a real file."""
    _, prodrg_param = get_prodrg_exec()

    if prodrg_param is not None:
        assert prodrg_param.is_file(), f"{prodrg_param} is not a file"


def test_prodrg_binary_is_executable_when_available():
    """If prodrg is available, the binary should be executable."""
    prodrg_exec, _ = get_prodrg_exec()

    if prodrg_exec is not None:
        assert os.access(prodrg_exec, os.X_OK), f"{prodrg_exec} is not executable"


def test_get_prodrg_exec_returns_consistent_pair():
    """Both return values should be either both None or both not None."""
    prodrg_exec, prodrg_param = get_prodrg_exec()

    assert (prodrg_exec is None) == (prodrg_param is None)


def test_get_prodrg_exec_returns_none_when_no_binary_no_env(monkeypatch):
    """Returns (None, None) when binary missing and PRODRG_EXEC not set."""
    import haddock.libs.libutil as libutil_module

    monkeypatch.setattr(libutil_module, "get_arch", lambda: "nonexistent-arch")
    monkeypatch.delenv("PRODRG_EXEC", raising=False)

    result = get_prodrg_exec()

    assert result == (None, None)


def test_get_prodrg_exec_uses_env_var(monkeypatch, tmp_path):
    """Uses PRODRG_EXEC env var when bundled binary is not found."""
    import haddock.libs.libutil as libutil_module

    fake_binary = tmp_path / "prodrg_custom"
    fake_binary.touch()
    monkeypatch.setattr(libutil_module, "get_arch", lambda: "nonexistent-arch")
    monkeypatch.setenv("PRODRG_EXEC", str(fake_binary))

    prodrg_exec, prodrg_param = get_prodrg_exec()

    assert prodrg_exec == fake_binary
    assert prodrg_param is not None
    assert prodrg_param.name == "prodrg.param"
