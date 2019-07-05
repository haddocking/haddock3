import itertools
import subprocess
import os

from haddock.workflows.scoring.pdb.utils import load_structure
from haddock.workflows.scoring.config import load_parameters

param_dic = load_parameters()
dcomplex_exe = param_dic['third-party']['dcomplex_exe']


# TODO: create Dfire symbolic links and correct locations
# def setup_dfire():

def dfire(pdb_f):

	prot_dic = load_structure(pdb_f)
	dfire_dic = {'binding': [], 'score': []}
	for segid_x, segid_y in itertools.combinations(prot_dic, 2):
		with open('dfire.inp','w') as inp:
			inp.write(f'{pdb_f}\n1\n{segid_x}\n1\n{segid_y}')
		inp.close()

		inp = open('dfire.inp')
		p = subprocess.run(dcomplex_exe, stdin=inp, stdout=subprocess.PIPE)
		try:
			_, binding, dfire_score = p.stdout.decode('utf-8').split()
			binding = float(binding.split('kcal/mol')[0])
			float(dfire_score)
		except ValueError:
			binding = float('nan')
			dfire_score = float('nan')

		dfire_dic['binding'].append(binding)
		dfire_dic['score'].append(dfire_score)

		os.remove('dfire.inp')

	binding_l = [float(e) for e in dfire_dic['binding']]
	score_l = [float(e) for e in dfire_dic['score']]

	mean_binding = sum(binding_l) / len(binding_l)
	mean_score = sum(score_l) / len(score_l)

	return mean_binding, mean_score
