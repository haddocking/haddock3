import logging
import io
import os
import sys
from functools import partial
from pathlib import Path


info_name = 'haddock3.log'
debug_name = 'haddock3.debug'

info_formatter = '[%(asctime)s]%(message)s'
debug_formatter = \
    "[%(asctime)s]%(filename)s:%(name)s:%(funcName)s:%(lineno)d: %(message)s"

log_levels = {
    'DEBUG': logging.DEBUG,
    'INFO': logging.INFO,
    'WARNING': logging.WARNING,
    'ERROR': logging.ERROR,
    'CRITICAL': logging.CRITICAL,
    }


def add_loglevel_arg(parser):
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=tuple(log_levels.keys()),
        )
    return


def add_streamhandler(
        log,
        stream,
        log_level='INFO',
        formatter=info_formatter,
        ):
    """Add a StreamHandler to the log object."""
    ch = logging.StreamHandler(stream)
    ch.setLevel(log_levels[log_level.upper()])
    ch.setFormatter(logging.Formatter(formatter))
    log.addHandler(ch)
    return ch


def set_log_for_cli(log, log_level, stream_to_stdout=True):
    """Set the logging streams for HADDOCK3 main client"""
    log.setLevel(log_levels[log_level])

    if log_level.upper() in ("INFO", "WARNING", "ERROR", "CRITICAL"):
        h = add_info_stringio(log)
        return [h.stream], [info_name]

    elif log_level.upper() in ("DEBUG"):
        if stream_to_stdout:
            log.handlers.clear()
            add_sysout_handler(log, log_level='DEBUG', formatter=debug_formatter)
        h1  = add_info_stringio(log)
        h2  = add_debug_stringio(log)
        return (h1.stream, h2.stream), (info_name, debug_name)

    else:
        raise AssertionError("Execution shouldn't reach this point.")


add_sysout_handler = partial(add_streamhandler, stream=sys.stdout)
add_stringio_handler = partial(add_streamhandler, stream=io.StringIO(newline=os.linesep))
add_info_stringio = partial(add_stringio_handler, log_level='INFO', formatter=info_formatter)
add_debug_stringio = partial(add_stringio_handler, log_level='DEBUG', formatter=debug_formatter)
