"""Handles configuration of the library"""
import json
import configparser
import sys
from utils.files import get_full_path
import os

# Locate folder and configuration files

etc_folder = get_full_path('haddock', 'etc')


config_file = os.path.join(etc_folder, 'scoring.ini')

environment = 'alcazar'
ini = configparser.ConfigParser(os.environ)
ini.read(config_file, encoding='utf-8')

scoring_json = None
if len(sys.argv) > 1:
    # blah = sys.argv[1]
    scoring_json = open(sys.argv[1])
else:
    print('ERROR: scoring.json not found\n')
    # exit()


try:
    # param_dic = json.loads(open(args.json_file).read())
    param_dic = json.loads(scoring_json.read())
except:
    print('ERROR: Check the json format\n')
    # exit()
