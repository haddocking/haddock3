"""Manage logging and logging helper functions."""
import io
import logging
import sys
from argparse import ArgumentParser
from functools import partial
from logging import FileHandler, Logger, StreamHandler
from typing import TextIO, Union

from haddock.core.typing import FilePath, StreamHandlerT


log_file_name = 'log'


info_formatter = '[%(asctime)s %(module)s %(levelname)s] %(message)s'
debug_formatter = (
    "[%(asctime)s] "
    "%(filename)s:%(name)s:%(funcName)s:%(lineno)d: "
    "%(message)s"
    )

log_formatters = {
    'DEBUG': debug_formatter,
    'INFO': info_formatter,
    'WARNING': info_formatter,
    'ERROR': info_formatter,
    'CRITICAL': info_formatter,
    }

log_levels = {
    'DEBUG': logging.DEBUG,
    'INFO': logging.INFO,
    'WARNING': logging.WARNING,
    'ERROR': logging.ERROR,
    'CRITICAL': logging.CRITICAL,
    }


def add_loglevel_arg(parser: ArgumentParser) -> None:
    """Add log level argument to CLI."""
    parser.add_argument(
        "--log-level",
        default='INFO',
        choices=list(log_levels.keys()),
        )
    return


def add_handler(
        log_obj: Logger,
        handler: type[StreamHandlerT],
        stream: Union[FilePath, TextIO],
        log_level: str = 'INFO',
        formatter: str = info_formatter,
        ) -> StreamHandlerT:
    """Add a logging Handler to the log object."""
    ch = handler(stream)
    ch.setLevel(log_levels[log_level.upper()])
    ch.setFormatter(logging.Formatter(formatter))
    log_obj.addHandler(ch)
    return ch


def add_log_for_CLI(log: Logger, log_level: str, logfile: FilePath) -> None:
    """Configure log for command-line clients."""
    llu = log_level.upper()

    params = {
        'log_level': llu,
        'formatter': log_formatters[llu],
        }

    log.handlers.clear()
    add_sysout_handler(log, **params)
    add_logfile_handler(log, stream=logfile, **params)
    return


add_sysout_handler = partial(add_handler, handler=StreamHandler, stream=sys.stdout)  # noqa: E501
add_syserr_handler = partial(add_handler, handler=StreamHandler, stream=sys.stderr, log_level='ERROR')  # noqa: E501
add_logfile_handler = partial(add_handler, handler=FileHandler, stream=log_file_name)  # noqa: E501
add_stringio_handler = partial(add_handler, handler=StreamHandler, stream=io.StringIO())  # noqa: E501
