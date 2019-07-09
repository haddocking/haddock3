# """Handles configuration of the library and the parameter loading"""
import json
import sys
import platform

py_version = float(platform.python_version()[:3])
if py_version < 3.6:
	raise SystemExit("ERROR: required Python version is 3.6.x")


def load_parameters():
	# TODO: check if all required parameters are filled in
	param_dic = None
	scoring_json = open(sys.argv[1])
	try:
		param_dic = json.loads(scoring_json.read())
	except:
		print('ERROR: scoring.json not found\n')
		exit()

	return param_dic
