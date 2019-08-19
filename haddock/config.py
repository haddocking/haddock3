"""Handles configuration of the library"""
import json
import configparser
import sys
from utils.files import get_full_path
import os

# Locate folder and configuration files

etc_folder = get_full_path('haddock', 'etc')
config_file = os.path.join(etc_folder, 'haddock3.ini')

# environment = 'alcazar'
ini = configparser.ConfigParser(os.environ)
ini.read(config_file, encoding='utf-8')

run_json = open(sys.argv[1])

try:
    # param_dic = json.loads(open(args.json_file).read())
    param_dic = json.loads(run_json.read())
except:
    print('ERROR: scoring.json not found\n')
    exit()


# """Handles configuration of the library"""
#
# import configparser
# import os
# from pathlib import PosixPath
# from haddock.modules.error import HaddockError
#
#
# # Locate folder and configuration files
# etc_folder = os.path.join(os.path.dirname(__file__), 'etc')
# config_file = os.path.join(etc_folder, 'haddock3.conf')
#
# parser = configparser.ConfigParser()
# parser.read(config_file, encoding='utf-8')
#
# # Intended for different environments in the future
# environment = 'development'
#
# # Preferred method is container
# CNS_bin = None
# host_data_folder = PosixPath(parser.get(environment, 'host_data_folder')).expanduser()
# container_data_folder = parser.get(environment, 'container_data_folder')
# CNS_container = parser.get(environment, 'CNS_container')
# if not CNS_container:
#     CNS_bin = parser.get(environment, 'CNS_bin')
#     if not CNS_bin:
#         try:
#             CNS_bin = os.environ['CNS_BIN']
#         except KeyError:
#             raise HaddockError('CNS container not found or system variable CNS_BIN is not defined')
# else:
#     CNS_container = CNS_container.format(host_data_folder, container_data_folder)

