import logging


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
    ch = logging.StreamHandler()
    ch.setLevel(log_levels[log_level.upper()])
    ch.setFormatter(logging.Formatter(formatter))
    log.addHandler(ch)
    return


def set_log_level_for_clis(log, log_level, add_streamhandler=True):
    """."""
    log.setLevel(log_levels[log_level])

    if log_level.upper() in ("INFO", "WARNING", "ERROR", "CRITICAL"):
        init_info_file(log)

    elif log_level.upper() in ("DEBUG"):
        log.handlers.clear()
        if add_streamhandler:
            add_streamhandler(log, log_level='DEBUG', formatter=debug_formatter)
        init_log_files(log)
        init_debug_file(log)

    return


def init_debug_file(
        log,
        file_name='haddock3.debug',
        formatter=debug_formatter,
        ):
    """."""
    handler = logging.FileHandler(file_name, mode='w')
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(formatter))
    log.addHandler(handler)
    return


def init_info_file(
        log,
        file_name='haddock3.log',
        formatter=info_formatter,
        ):
    """
    """
    handler = logging.FileHandler(file_name, mode='w')
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(formatter))
    log.addHandler(handler)
    return
