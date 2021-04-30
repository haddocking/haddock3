"""Command line messages"""
import sys
import os
import random
import logging
from datetime import datetime
from haddock.version import CURRENT_VERSION

logger = logging.getLogger(__name__)


def greeting():
    """Initial message"""
    now = datetime.now().replace(second=0, microsecond=0)
    python_version = sys.version
    message = (f"""{os.linesep}##############################################{os.linesep}"""
               f"""#                                            #{os.linesep}"""
               f"""#                  HADDOCK3                  #{os.linesep}"""
               f"""#                                            #{os.linesep}"""
               f"""##############################################{os.linesep}"""
               f"""{os.linesep}"""
               f"""Starting HADDOCK {CURRENT_VERSION} on {now}{os.linesep}"""
               f"""{os.linesep}"""
               f"""Python {python_version}{os.linesep}"""
               )
    return message


def get_greetings(how_many=3):
    """Get different (how_many) greeting messages"""
    greetings = ["Tot ziens!", "Good bye!", "Até logo!", "Ciao!", "Au revoir!", "Adéu-siau!", "Agur!", "Dovidenia!"]
    n = how_many % len(greetings)
    return " ".join(random.sample(greetings, k=n))


def adieu():
    """Final message"""
    end = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    bye = get_greetings()
    message = (f"""Finished at {end}. {bye}""")
    return message
