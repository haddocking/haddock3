# """Handles configuration of the library and the parameter loading"""

# import configparser
import json
import sys
import os
import platform

# from pathlib import PosixPath
# from haddock.modules.error import HaddockError


py_version = float(platform.python_version()[:3])
if py_version < 3.6:
	raise SystemExit("ERROR: required Python version is 3.6.x")


def load_parameters():
	# TODO: check if all required parameters are filled in
	param_dic = None
	# scoring_json = open('scoring.json')
	scoring_json = open(sys.argv[1])
	try:
		param_dic = json.loads(scoring_json.read())
	except:
		print('ERROR: scoring.json not found\n')
		exit()

	# param_dic['wd'] = os.path.dirname(__file__)
	return param_dic

# # Locate folder and configuration files
# etc_folder = os.path.join(os.path.dirname(__file__), 'etc')
# config_file = os.path.join(etc_folder, 'scoring.conf')
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
# # CNS_container = parser.get(environment, 'CNS_container')
# # if not CNS_container:
# CNS_bin = parser.get(environment, 'CNS_bin')
# 	if not CNS_bin:
# 		try:
# 			CNS_bin = os.environ['CNS_BIN']
# 		except KeyError:
# 			print('Oh no!')
# 		# raise HaddockError('CNS container not found or system variable CNS_BIN is not defined')
# else:
# 	CNS_container = CNS_container.format(host_data_folder, container_data_folder)
