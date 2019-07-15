"""Handles configuration of the library"""
import json
import os
import configparser
import argparse
# from pathlib import PosixPath
# from haddock.modules.error import HaddockError

# Locate folder and configuration files
etc_folder = os.path.join(os.path.dirname(__file__), 'etc')
config_file = os.path.join(etc_folder, 'scoring.ini')

environment = 'alcazar'
ini = configparser.ConfigParser()
ini.read(config_file, encoding='utf-8')

parser = argparse.ArgumentParser()
parser.add_argument("json_file", help="your json file")
args = parser.parse_args()

# scoring_json = open(sys.argv[1])
try:
    param_dic = json.loads(open(args.json_file).read())
except:
    print('ERROR: scoring.json not found\n')
    exit()


# print(args.echo)

# Intended for different environments in the future
# environment = 'alcazar'

# from configparser import ConfigParser
# import os
#
#
# class Config:
#     """Interact with configuration variables."""
#
#     configParser = ConfigParser()
#     etc_folder = os.path.join(os.path.dirname(__file__), 'etc')
#     config_file = os.path.join(etc_folder, 'scoring.ini')
#     # configFilePath = (os.path.join(os.getcwd(), 'scoring.ini'))
#
#     @classmethod
#     def initialize(cls, newhire_table):
#         """Start config by reading config.ini."""
#         cls.configParser.read(cls.config_file)
#
#     @classmethod
#     def prod(cls, key):
#         """Get prod values from config.ini."""
#         return cls.configParser.get('PROD', key)
#
#     @classmethod
#     def dev(cls, key):
#         """Get dev values from config.ini."""
#         return cls.configParser.get('DEV', key)