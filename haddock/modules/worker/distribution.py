import os
# import config as config

string = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'


class JobCreator:

	def __init__(self):
		self.dic = {}

		# self.wd = os.getcwd()
		self.job_dir = 'jobs'
		self.out_dir = 'output'
		# if not os.path.isdir(f"{self.wd}/run{config.param_dic['input']['run']}jobs"):
		# 	os.system(f'mkdir {self.wd}/jobs')

		# if not os.path.isdir(f'{self.wd}/out'):
		# 	os.system(f'mkdir {self.wd}/out')

	def delegate(self, job_num, job_id, recipe_str, input_model):
		""" Give one recipe to input Model """

		# recipe_str = input_script.recipe
		# for i, model in enumerate(input_model_list):
		# model_name = model.split('/')[-1].split('.')[0]
		model_root = input_model.split('/')[-1].split('.')[0]
		job_str = '0' * (2 - len(str(job_num))) + str(job_num)
		input_f = f'{self.job_dir}/{job_id}_{job_str}.inp'
		output_f = f'{self.out_dir}/{job_id}_{job_str}.out'

		# if not os.path.isfile(output_f):

		tbw = '! Input structure\n'
		tbw += f'eval ($file="{input_model}")\n'
		tbw += f'eval ($file_root="{model_root}")\n'
		tbw += recipe_str

		with open(input_f, 'w') as f:
			f.write(tbw)
		f.close()

		self.dic[job_num] = (input_f, output_f)

		return self.dic

	def delegate_it0(self, job_num, job_id, recipe_str, input_psf, input_pdb):
		""" Give one recipe to this input model list """

		# This is exclusive for it0
		job_str = '0' * (6 - len(str(job_num))) + str(job_num)
		input_f = f'{self.job_dir}/{job_id}_{job_str}.inp'
		output_f = f'{self.out_dir}/{job_id}_{job_str}.out'

		tbw = '! Input Topologies\n'
		tbw += 'structure\n'
		for psf in input_psf:
			tbw += f'  @@{psf}\n'
		tbw += 'end\n'

		tbw += '\n'
		tbw += '! Input Molecules\n'
		for pdb in input_pdb:
			tbw += f'coor @@{pdb}\n'

		ncomp = len(input_pdb)
		tbw += f'eval ($ncomponents={ncomp})\n'

		for i in range(ncomp):
			tbw += f'eval ($prot_segid_mol{i+1}="{string[i]}")\n'

		# model_root = input_pdb[0].split('/')[-1].split('.')[0]
		# tbw += f'eval ($file_root="{model_root}")\n'

		tbw += recipe_str

		with open(input_f, 'w') as f:
			f.write(tbw)
		f.close()

		self.dic[job_num] = (input_f, output_f)

		return self.dic
