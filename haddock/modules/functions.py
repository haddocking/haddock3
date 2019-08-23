# placeholder
import os
import config as config


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
		print(f'ERROR: {run_dir} already present')
		# exit()

	if not os.path.isdir(data_dir):
		os.system(f'mkdir {data_dir}')
	for mol_id in config.param_dic['input']['molecules']:
		molecule = config.param_dic['input']['molecules'][mol_id]

		if molecule == mol_id + '.pdb':
			print("ERROR: Name mol_X.pdb not supported, please rename your molecule.")
			exit()

		os.system(f'cp {molecule} {run_dir}/data/{mol_id}.pdb')
		config.param_dic['input']['molecules'][mol_id] = f'data/{mol_id}.pdb'

	if not os.path.isdir(begin_dir):
		os.system(f'mkdir {begin_dir}')

	if not os.path.isdir(job_dir):
		os.system(f'mkdir {job_dir}')

	if not os.path.isdir(output_dir):
		os.system(f'mkdir {output_dir}')

	os.chdir(run_dir)


def retrieve_output(jobs):
	# First attempt, parse .out and identify the tag
	output_dic = {}
	for j in jobs.dic:
		_, out = jobs.dic[j]
		output_dic[j] = []
		for line in reversed(list(open(out))):
			if 'OUTPUT:' in line and not 'CNS' in line:
				v = line.rstrip().split(':')[-1]
				output_dic[j].append(v)
	return output_dic
