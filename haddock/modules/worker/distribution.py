import os


class JobCreator:

	def __init__(self, job_id, job_folder):
		self.dic = {}
		self.folder = job_folder
		self.id = job_id

	def delegate(self, job_num, input_file_str):

		if self.folder:
			path = self.folder
		else:
			path = os.getcwd()

		if not os.path.isdir(path):
			os.mkdir(path)

		job_str = '0' * (6 - len(str(job_num))) + str(job_num)

		input_f = f'{path}/{self.id}_{job_str}.inp'
		output_f = f'{path}/{self.id}_{job_str}.out'

		# Q: Would it be a good idea to change this to a paralel implementation? or would it stress the I/O too much?
		#  For example for it0 the input file is ~11000 lines and ~350K
		with open(input_f, 'w') as f:
			f.write(input_file_str)
		f.close()

		self.dic[job_num] = (input_f, output_f)

		return self.dic
