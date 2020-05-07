import itertools
import subprocess
import os
import re
from haddock.utils.files import get_full_path
from haddock.modules.structure.utils import PDB


def dfire(pdb_f, dcomplex_exe):

	fixed_pdb = PDB.fix_id(pdb_f)
	prot_dic = PDB.load_structure(fixed_pdb)

	dfire_regex = r'-?\d*?\.\d\d*'

	dfire_dic = {'binding': [], 'score': []}
	for segid_x, segid_y in itertools.combinations(prot_dic, 2):

		input_str = str.encode(f'{fixed_pdb}\n1\n{segid_x}\n1\n{segid_y}')

		p = subprocess.run(dcomplex_exe, input=input_str, stdout=subprocess.PIPE)

		try:
			result = p.stdout.decode('utf-8')
			score, binding = re.findall(dfire_regex, result)

			score = float(score)
			binding = float(binding)

		except ValueError:

			binding = .0
			score = .0

		dfire_dic['binding'].append(binding)
		dfire_dic['score'].append(score)

	binding_l = [float(e) for e in dfire_dic['binding']]
	score_l = [float(e) for e in dfire_dic['score']]

	mean_binding = sum(binding_l) / len(binding_l)
	mean_score = sum(score_l) / len(score_l)

	return mean_binding, mean_score
