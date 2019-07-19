import multiprocessing
import subprocess
import haddock.workflows.scoring.config as config


class CNS:
	"""Wraps interacting with the CNS core facilities"""

	def __init__(self):
		self.stdout = None
		self.cns_exec = config.ini.get('cns', 'cns_exe')
		self.nproc = config.param_dic['input']['nproc']
		self.run_scheme = config.param_dic['input']['run_scheme']

	def run(self, jobs):
		"""Pass the input defined in input_file to CNS"""
		if self.run_scheme == 'serial':
			self.run_serial(jobs)

		elif self.run_scheme == 'alcazar':
			self.run_alcazar(jobs)

		elif self.run_scheme == 'parallel':
			self.run_parallel(jobs)

		else:
			print('ERROR: run_scheme not recognized')
			exit()

	def run_serial(self, jobs):

		for job_id in jobs.job_dic:
			input_f = jobs.job_dic[job_id][0]
			output_f = jobs.job_dic[job_id][1]

			with open(input_f) as inp:
				with open(output_f, 'w+') as outf:

					p = subprocess.Popen(self.cns_exec, stdin=inp, stdout=outf, close_fds=True)

					out, error = p.communicate()

					p.kill()

					if error:
						print('Oh no!')
						exit()
					# raise CNSError(error)

		return out

	def run_alcazar(self, jobs):
		pass

	def run_parallel(self, jobs, nproc=48):
		""" Execute the jobs in parallel """
		print(f'+ Running parallel, {nproc} cores')
		arg_list = []
		for job_id in jobs.job_dic:
			inp_f = jobs.job_dic[job_id][0]
			out_f = jobs.job_dic[job_id][1]
			arg_list.append((self.cns_exec, inp_f, out_f))

		with multiprocessing.Pool(processes=nproc) as pool:
			_ = pool.starmap(self.execute, arg_list)

	@staticmethod
	def execute(cns_exec, input_f, output_f):
		with open(input_f) as inp:
			with open(output_f, 'w+') as outf:
				p = subprocess.Popen(cns_exec, stdin=inp, stdout=outf, close_fds=True)
				_, _ = p.communicate()
				p.kill()





