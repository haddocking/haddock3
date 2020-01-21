import configparser
import subprocess
import os
from haddock.utils.files import get_full_path

etc_folder = get_full_path('haddock', 'etc')
config_file = os.path.join(etc_folder, 'haddock3.ini')

ini = configparser.ConfigParser(os.environ)
ini.read(config_file, encoding='utf-8')

reduce_exec = ini.get('third party', 'reduce_exe')


class ReduceError(Exception):
	"""
    Custom exception class to provide better error messages when calling reduce from molprobity.
    """

	def __init__(self, user_message):
		# Add information to our message

		full_message = "There was an error while calling Reduce.\n{}".format(user_message)
		# Call base class constructor
		Exception.__init__(self, full_message)


def run_reduce(pdb_f):
	""""
	Reads a PDB file and outputs the corrected structure and a dictionary with protonation states.
	Expects either an open file handle or a string with a PDB formatted structure.

	Option strip_header removes all lines not starting with ATOM, TER, END, etc.. (check PDB format)
	"""

	cmd = [reduce_exec, '-BUILD', '-Xplor', '-quiet', pdb_f]

	try:
		process_handle = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
		                                  stderr=subprocess.PIPE, close_fds=False)
	except Exception as e:
		raise ReduceError("There was an error running the 'reduce': %s\n(%s)\n" % (' '.join(cmd), e))

	p_stdout, p_stderr = process_handle.communicate()

	return p_stdout, p_stderr


def analyze_protonation_state(pdb_f):

	out, err = run_reduce(pdb_f)
	out = out.decode('utf-8')

	his_db = {}
	hisprotonatoms = {" HD1": 'd1', " HD2": 'd2', " HE1": 'e1', " HE2": 'e2'}

	# Build Histidine 'Database'
	# Assign protonation states based on presence of atoms in the PDB file
	for l in out.split('\n'):
		if l.startswith('ATOM'):
			chain = l[21]
			aname = l[12:16]
			resn = l[17:20]
			resi = int(l[22:26])

			try:
				_ = his_db[chain]
			except:
				his_db[chain] = {}

			if resn != "HIS":
				continue
			elif resn == "HIS" and aname in hisprotonatoms:
				if resi not in his_db[chain]:
					his_db[chain][resi] = {}
				currhis = his_db[chain][resi]
				histidine_state = hisprotonatoms.get(aname)
				currhis[histidine_state] = True
			else:
				if resi not in his_db[chain]:
					his_db[chain][resi] = {}
	# Decide on Protonation State for CNS/HADDOCK
	ret = {}
	for chain in his_db:
		ret[chain] = {}
		for resi in his_db[chain]:
			if his_db[chain][resi]:
				his = his_db[chain][resi]
				dcount = his.get('d1', 0) + his.get('d2', 0)
				ecount = his.get('e1', 0) + his.get('e2', 0)
				total_count = dcount + ecount

				# To easily retrive histidines protonation states we take resi as key and the state as value [{132: "HISE"}]
				if total_count == 4:
					ret[chain][resi] = "HIS+"
				elif total_count == 3:
					if dcount == 2:
						ret[chain][resi] = "HISD"
					else:
						ret[chain][resi] = "HISE"
				else:
					raise ReduceError("MolProbity could not guess the protonation state of histidine %d in %s"
					                      % (resi, pdb_f))
			else:
				ret[chain][resi] = None
	return ret
