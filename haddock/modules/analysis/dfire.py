import configparser
import itertools
import subprocess
import os
from haddock.utils.files import get_full_path
from haddock.modules.structure.utils import PDB

etc_folder = get_full_path('haddock', 'etc')
config_file = os.path.join(etc_folder, 'haddock3.ini')

ini = configparser.ConfigParser(os.environ)
ini.read(config_file, encoding='utf-8')

dcomplex_exe = ini.get('third_party', 'dcomplex_exe')


def dfire(pdb_f):

	prot_dic = PDB.load_structure(pdb_f)
	dfire_dic = {'binding': [], 'score': []}
	for segid_x, segid_y in itertools.combinations(prot_dic, 2):

		p = subprocess.run([dcomplex_exe, pdb_f, segid_x, segid_y], stdout=subprocess.PIPE)

		try:
			binding = float(p.stdout.decode('utf-8').split()[-2])
		except ValueError:
			binding = float('nan')

		dfire_dic['binding'].append(binding)

	binding_l = [float(e) for e in dfire_dic['binding']]
	mean_binding = sum(binding_l) / len(binding_l)

	return mean_binding
