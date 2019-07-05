import os


class JobCreator:

	def __init__(self):
		self.job_dic = {}
		self.wd = os.getcwd()

		if not os.path.isdir(f'{self.wd}/jobs'):
			os.system(f'mkdir {self.wd}/jobs')

		if not os.path.isdir(f'{self.wd}/out'):
			os.system(f'mkdir {self.wd}/out')

	def delegate(self, input_script, input_model_list):
		""" Give one recipe to each Model """

		recipe_str = input_script.recipe
		for i, model in enumerate(input_model_list):
			model_name = model.split('/')[1].split('.')[0]

			input_f = f'{self.wd}/jobs/{model_name}.inp'
			output_f = f'{self.wd}/out/{model_name}.out'

			tbw = f'eval($file="{model}")\n' + recipe_str
			with open(input_f, 'w') as f:
				f.write(tbw)
			f.close()

			self.job_dic[i] = (input_f, output_f)
