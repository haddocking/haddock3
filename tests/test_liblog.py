"""Test liblog."""

import argparse
import tempfile

import pytest

from haddock import log
from haddock.libs import liblog


def test_add_systout_handler():
    log.info("something")
    assert len(log.handlers) == 1


def test_add_logfile_handler(monkeypatch):
    with tempfile.NamedTemporaryFile(delete=False) as temp_f:
        log_fname = temp_f.name
    monkeypatch.setattr("haddock.libs.liblog.log_file_name", log_fname)
    log.info("something")
    assert len(log.handlers) == 1


def test_dics_keys():
    """Test log levels are okay in dics."""
    a = list(liblog.log_formatters.keys())
    b = list(liblog.log_levels.keys())

    assert a == b
    assert len(a) == 5


@pytest.mark.parametrize("key", liblog.log_levels.keys())
def test_add_loglevel_arg(key):
    """Test adds log level argument to CLI."""
    ap = argparse.ArgumentParser()
    liblog.add_loglevel_arg(ap)
    cmd = ap.parse_args(f"--log-level {key}".split())
    assert cmd.log_level


def test_add_loglevel_arg_error():
    """Test adds log level argument to CLI."""
    ap = argparse.ArgumentParser()
    liblog.add_loglevel_arg(ap)
    with pytest.raises(SystemExit) as err:
        ap.parse_args("--log-level BAD".split())

    assert err.value.code == 2
