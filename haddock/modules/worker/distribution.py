import os


class JobCreator:

	def __init__(self):
		self.dic = {}

	def delegate(self, job_num, job_id, recipe_str, job_folder=None):

		if job_folder:
			path = job_folder
		else:
			path = os.getcwd()

		if not os.path.isdir(path):
			os.mkdir(path)

		job_str = '0' * (6 - len(str(job_num))) + str(job_num)

		input_f = f'{path}/{job_id}_{job_str}.inp'
		output_f = f'{path}/{job_id}_{job_str}.out'

		# Would it be a good idea to change this to a paralel implementation? or would it stress the I/O too much?
		#  For example for it0 the input file is ~11000 lines and ~350K
		with open(input_f, 'w') as f:
			f.write(recipe_str)
		f.close()

		self.dic[job_num] = (input_f, output_f)

		return self.dic
