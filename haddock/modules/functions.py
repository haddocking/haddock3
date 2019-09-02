# placeholder
import os
import config as config

from haddock.modules.structure.utils import PDB


def init():
	# move input molecules to correct places
	run_dir = f"run{config.param_dic['input']['run']}"
	data_dir = f'{run_dir}/data'
	begin_dir = f'{run_dir}/begin'

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

	try:
		ambig_fname = config.param_dic['input']['restraints']['ambig']
		os.system(f'cp {ambig_fname} {run_dir}/data/ambig.tbl')
		config.param_dic['input']['restraints']['ambig'] = 'data/ambig.tbl'
	except KeyError:
		pass

	try:
		unambig_fname = config.param_dic['input']['restraints']['unambig']
		os.system(f'cp {unambig_fname} {run_dir}/data/unambig.tbl')
		config.param_dic['input']['restraints']['unambig'] = 'data/unambig.tbl'
	except KeyError:
		pass

	try:
		hbond_fname = config.param_dic['input']['restraints']['hbond']
		os.system(f'cp {hbond_fname} {run_dir}/data/hbond.tbl')
		config.param_dic['input']['restraints']['hbond'] = 'data/hbond.tbl'
	except KeyError:
		pass
	try:
		dihe_fname = config.param_dic['input']['restraints']['dihedrals']
		os.system(f'cp {dihe_fname} {run_dir}/data/dihe.tbl')
		config.param_dic['input']['restraints']['dihe_fname'] = 'data/dihe.tbl'
	except KeyError:
		pass

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
	# FIXME: This must read the output file and retrieve OUTPUT: and identify if the job was successful or not
	output_dic = {}
	complete_check = None
	for j in jobs.dic:
		inp, out = jobs.dic[j]
		output_dic[j] = []
		for i, line in enumerate(reversed(list(open(out)))):
			if 'OUTPUT:' in line and 'CNS' not in line:
				v = line.strip().split(':')[1][1:]
				output_dic[j].append(v)
			if 'CNSsolve>stop' in line:
				# passed ok
				complete_check = True
			if not complete_check and 'ABORT' in line:
				# job has failed!
				print(f'+ ERROR: Job {inp} has failed, please check the output {out}')
				exit()
			if i == 50:
				break
	return output_dic
