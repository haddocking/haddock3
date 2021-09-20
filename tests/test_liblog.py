"""Test liblog"""
import io
import logging

import pytest


from haddock import log
from haddock.libs import liblog


@pytest.mark.parametrize(
    'func',
    [
        liblog.add_sysout_handler,
        liblog.add_stringio_handler,
        liblog.add_info_stringio,
        liblog.add_debug_stringio,
        ]
    )
def test_add_sysout_handler(func):
    """."""
    log.handlers.clear()
    func(log)
    assert len(log.handlers) == 1


def test_set_log_for_cli():
    s, n = liblog.set_log_for_cli(log, log_level='info')
    assert isinstance(s, list)
    assert isinstance(n, list)
    assert isinstance(s[0], io.StringIO)
    assert n[0] == liblog.info_name

