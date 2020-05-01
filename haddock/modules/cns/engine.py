import multiprocessing
import random
import subprocess
import os
import toml
from haddock.modules.functions import load_ini

ini = load_ini('haddock3.ini')


class CNS:
	"""Wraps interacting with the CNS core facilities"""

	def __init__(self):
		self.stdout = None
		self.cns_exec = ini.get('DEFAULT', 'cns_exe')

		self.run_param = toml.load('data/run.toml')
		self.run_scheme = self.run_param['execution_parameters']['scheme']

	def validate_cns(self):

		cns_script_child = subprocess.Popen([self.cns_exec], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		stdout, stderr = cns_script_child.communicate(b"display Valid\nstop")
		if 'Valid' in stdout.decode('utf-8'):
			return True
		else:
			raise False


	@staticmethod
	def replace_seed(jobs):
		for j in jobs.dic:
			inp_f, out = jobs.dic[j]
			print(inp_f)
			new_inp_f = inp_f.replace('.inp', '.inp_')
			with open(new_inp_f, 'w') as n_fh:
				with open(inp_f, 'r') as fh:
					for line in fh.readlines():
						if 'eval ($seed=' in line:
							seed = random.randint(100, 999)
							line = f'eval ($seed={seed})\n'
						n_fh.write(line)
			# overwrite
			os.rename(new_inp_f, inp_f)
			# remove .out
			os.remove(out)
		return jobs

	def run(self, jobs, retry=False):
		"""Pass the input defined in input_file to CNS"""
		if retry:
			jobs = self.replace_seed(jobs)

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

		for job_id in jobs.dic:
			input_f = jobs.dic[job_id][0]
			output_f = jobs.dic[job_id][1]

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

	def run_parallel(self, jobs):
		""" Execute the jobs in parallel """
		job_n = len(jobs.dic)
		nproc = self.run_param['execution_parameters']['nproc']
		if not nproc:
			print('+ ERROR: With run_scheme parallel you must define nproc')
			exit()
		if nproc > job_n:
			nproc = job_n
		if nproc > multiprocessing.cpu_count():
			nproc = multiprocessing.cpu_count()

		print(f'+ Running parallel, {nproc} cores {job_n} jobs')
		arg_list = []
		for job_id in jobs.dic:
			inp_f = jobs.dic[job_id][0]
			out_f = jobs.dic[job_id][1]
			if not os.path.isfile(out_f):
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

class CNSError(Exception):
	def __init__(self, *args):
		if args:
			self.message = args[0]
		else:
			self.message = None
		raise SystemExit('{0} '.format(self.message))
