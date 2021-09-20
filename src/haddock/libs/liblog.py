import logging
import io
import os
import sys
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
        log_level='INFO',
        formatter=info_formatter,
        ):
    """Add streamhandler to the log object."""
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_levels[log_level.upper()])
    ch.setFormatter(logging.Formatter(formatter))
    log.addHandler(ch)
    return


def set_log_level_for_clis(log, log_level, stream_to_stdout=True):
    """."""
    log.setLevel(log_levels[log_level])

    if log_level.upper() in ("INFO", "WARNING", "ERROR", "CRITICAL"):
        h, n = init_info_file(log)
        return [h.stream], [n]

    elif log_level.upper() in ("DEBUG"):
        log.handlers.clear()
        if stream_to_stdout:
            add_streamhandler(log, log_level='DEBUG', formatter=debug_formatter)
        h1, n1 = init_info_file(log)
        h2, n2 = init_debug_file(log)
        return (h1.stream, h2.stream), (n1, n2)

    else:
        raise AssertionError("Execution shouldn't reach this point.")


def init_debug_file(
        log,
        file_name=debug_name,
        formatter=debug_formatter,
        ):
    """."""
    #handler = logging.FileHandler(file_name, mode='w')
    handler = logging.StreamHandler(io.StringIO(newline=os.linesep))
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(formatter))
    log.addHandler(handler)
    return handler, file_name


def init_info_file(
        log,
        file_name=info_name,
        formatter=info_formatter,
        ):
    """
    """
    handler = logging.StreamHandler(io.StringIO(newline=os.linesep))
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(formatter))
    log.addHandler(handler)
    return handler, file_name
