"""Handles configuration of the library"""
import toml
import sys
from haddock.modules.functions import load_ini

ini = load_ini('haddock3.ini')

toml_f = sys.argv[1]

if not sys.argv[1]:
    print('+ERROR: Setup file not found')
else:
    try:
        param_dic = toml.load(toml_f)
    except:
        pass
        # print('+ERROR: There is something wrong with your setup file')
