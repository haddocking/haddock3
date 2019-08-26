# placeholder
import os
import config as config

from haddock.modules.structure.utils import PDB


def init():
	# move input molecules to correct places
	wd = os.getcwd()
	run_dir = f"run{config.param_dic['input']['run']}"
	data_dir = f'{run_dir}/data'
	begin_dir = f'{run_dir}/begin'
	job_dir = f'{run_dir}/jobs'
	output_dir = f'{run_dir}/output'

	if not os.path.isdir(run_dir):
		os.system(f'mkdir {run_dir}')
	else:
		print(f'+ WARNING: {run_dir} already present')
	# exit()

	if not os.path.isdir(data_dir):
		os.system(f'mkdir {data_dir}')
	for mol_id in config.param_dic['input']['molecules']:
		molecule = config.param_dic['input']['molecules'][mol_id]

		if molecule == mol_id + '.pdb':
			print("+ ERROR: Name mol_X.pdb not supported, please rename your molecule.")
			exit()

		# ensembles will be treated later
		os.system(f'cp {molecule} {run_dir}/data/{mol_id}_1.pdb')

		config.param_dic['input']['molecules'][mol_id] = f'data/{mol_id}_1.pdb'

	if not os.path.isdir(begin_dir):
		os.system(f'mkdir {begin_dir}')

	if not os.path.isdir(job_dir):
		os.system(f'mkdir {job_dir}')

	if not os.path.isdir(output_dir):
		os.system(f'mkdir {output_dir}')

	os.chdir(run_dir)


def treat_ensemble(pdb_dic):
	check = False
	new_d = {}
	for mol in pdb_dic:
		new_d[mol] = []
		pdb = pdb_dic[mol]
		with open(pdb) as f:
			reversed_list = reversed(f.readlines())
			for line in reversed_list:
				if 'ENDMDL' in line:
					check = True
					break
			if check:
				split_models = PDB.split_models(pdb)
				for model in split_models:
					new_d[mol].append(model)
				continue
			else:
				new_d[mol].append(pdb)

		f.close()
	return new_d


def molecularize(dic):
	out_dic = {}
	for e in dic:
		root = dic[e][0].split('/')[1].split('_')[0]
		try:
			_ = out_dic[root]
		except:
			out_dic[root] = []

		out_dic[root].append(tuple(dic[e]))

	return out_dic


def retrieve_output(jobs):
	# First attempt, parse .out and identify the tag
	output_dic = {}
	for j in jobs.dic:
		_, out = jobs.dic[j]
		output_dic[j] = []
		for line in reversed(list(open(out))):
			if 'OUTPUT:' in line and not 'CNS' in line:
				v = line.strip().split(':')[1][1:]
				output_dic[j].append(v)
	return output_dic
