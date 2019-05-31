"""Handles configuration of the library"""

import configparser
import os
from haddock.modules.error import HaddockError


# Locate folder and configuration files
etc_folder = os.path.join(os.path.dirname(__file__), 'etc')
config_file = os.path.join(etc_folder, 'haddock3.conf')

parser = configparser.ConfigParser()
parser.read(config_file, encoding='utf-8')

CNS_bin = parser.get('development', 'CNS_bin')
if not CNS_bin:
    try:
        CNS_bin = os.environ['CNS_BIN']
    except KeyError:
        raise HaddockError('CNS not found or system variable CNS_BIN is not defined')
