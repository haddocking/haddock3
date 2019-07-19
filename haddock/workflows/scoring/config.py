"""Handles configuration of the library"""
import json
import os
import configparser
import argparse
import sys
# from pathlib import PosixPath
# from haddock.modules.error import HaddockError

# Locate folder and configuration files
etc_folder = os.path.join(os.path.dirname(__file__), 'etc')
config_file = os.path.join(etc_folder, 'scoring.ini')

environment = 'alcazar'
ini = configparser.ConfigParser()
ini.read(config_file, encoding='utf-8')

# parser = argparse.ArgumentParser()
# parser.add_argument("json_file", help="your json file")
# args = parser.parse_args()

scoring_json = open(sys.argv[1])

try:
    # param_dic = json.loads(open(args.json_file).read())
    param_dic = json.loads(scoring_json.read())
except:
    print('ERROR: scoring.json not found\n')
    exit()
