"""Test liblog."""
import argparse
from logging import FileHandler, StreamHandler

import pytest

from haddock import log
from haddock.libs import liblog


@pytest.mark.parametrize(
    'func',
    [
        liblog.add_sysout_handler,
        liblog.add_logfile_handler,
        ],
    )
def test_add_handlers(func):
    """Test add handlers functions."""
    log.handlers.clear()
    func(log)
    log.info('something')
    assert len(log.handlers) == 1


def test_dics_keys():
    """Test log levels are okay in dics."""
    a = list(liblog.log_formatters.keys())
    b = list(liblog.log_levels.keys())

    assert a == b
    assert len(a) == 5


@pytest.mark.parametrize('key', liblog.log_levels.keys())
def test_add_loglevel_arg(key):
    """Test adds log level argument to CLI."""
    ap = argparse.ArgumentParser()
    liblog.add_loglevel_arg(ap)
    cmd = ap.parse_args(f'--log-level {key}'.split())
    assert cmd.log_level


def test_add_loglevel_arg_error():
    """Test adds log level argument to CLI."""
    ap = argparse.ArgumentParser()
    liblog.add_loglevel_arg(ap)
    with pytest.raises(SystemExit) as err:
        ap.parse_args('--log-level BAD'.split())

    assert err.value.code == 2


def test_close_logfile_handlers(tmp_path):
    """File handlers are closed and detached, stream handlers are kept."""
    log.handlers.clear()
    logfile = tmp_path / "log"
    liblog.add_sysout_handler(log)
    liblog.add_logfile_handler(log, stream=str(logfile))
    log.info("something")

    file_handlers = [h for h in log.handlers if isinstance(h, FileHandler)]
    assert len(file_handlers) == 1

    liblog.close_logfile_handlers(log)

    # No FileHandler remains, the StreamHandler is untouched.
    assert not any(isinstance(h, FileHandler) for h in log.handlers)
    assert any(isinstance(h, StreamHandler) for h in log.handlers)
    # The underlying file descriptor has been released.
    assert file_handlers[0].stream is None or file_handlers[0].stream.closed
