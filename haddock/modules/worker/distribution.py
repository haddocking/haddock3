# import os
# import config as config


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

		input_f = f'{self.job_dir}/{job_id}_{job_num}.inp'
		output_f = f'{self.out_dir}/{job_id}_{job_num}.out'

		tbw = '! Input structure\n'
		tbw += f'eval ($file="{input_model}")\n'
		tbw += f'eval ($file_root="{model_root}")\n'
		tbw += recipe_str

		with open(input_f, 'w') as f:
			f.write(tbw)
		f.close()

		self.dic[job_num] = (input_f, output_f)

		return self.dic
