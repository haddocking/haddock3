import os
import random
import glob
from haddock.modules.functions import load_ini
from haddock.modules.structure.utils import PDB

ini = load_ini('haddock3.ini')


class InputGenerator:

	def __init__(self, recipe_file, input_folder=''):

		with open(recipe_file, 'r') as fh:
			self.recipe_str = ''.join(fh.readlines())
		fh.close()

		self.folder = input_folder + '/'

	def generate(self, protonation_dic, output_pdb, output_psf, input_pdb, input_psf, output_fname):

		forcefield_parameters = ini.get('parameters', 'param_file')
		param = self.load_ff_parameters(forcefield_parameters)

		forcefield_topology = ini.get('topology', 'top_file')
		top = self.load_ff_topology(forcefield_topology)

		molecule_linkage = ini.get('link', 'link')
		link = self.load_link(molecule_linkage)

		protonation = self.load_protonation_state(protonation_dic)

		translation_vectors = dict([(i, ini.get('translation_vectors', f'trans_vector_{i}')) for i in range(51)])
		trans_vec = self.load_trans_vectors(translation_vectors)

		tensor_options = ['tensor_psf', 'tensor_pdb', 'tensor_para_psf', 'tensor_para_pdb', 'tensor_dani_psf']
		tensor_params = dict([(e, ini.get('tensor', e)) for e in tensor_options])
		tensor = self.load_tensor(tensor_params)

		scatter_lib = ini.get('scatter', 'scatter_lib')
		scatter = self.load_scatter(scatter_lib)

		axis_options = ['top_axis', 'par_axis', 'top_axis_dani']
		axis_params = dict([(e, ini.get('axis', e)) for e in axis_options ])
		axis = self.load_axis(axis_params)

		output = self.prepare_output(output_pdb, output_psf, output_fname)
		input_str = self.prepare_input(input_pdb, input_psf)

		inp = param + top + input_str + output + link + protonation + trans_vec + tensor + scatter + axis + self.recipe_str

		return inp

	@staticmethod
	def prepare_input(pdb_input, psf_input=None):
		""" Write input of recipe """
		# This section will be written for any recipe
		#  Even if some CNS variables are not used, it should not be an issue.

		input_str = '\n! Input structure\n'

		string = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		ncomp = None

		if psf_input:
			if type(psf_input) == str:
				input_str += 'structure\n'
				input_str += f'  @@{psf_input}\n'
				input_str += 'end\n'
			if type(psf_input) == list:
				input_str += 'structure\n'
				for psf in psf_input:
					input_str += f'  @@{psf}\n'
				input_str += 'end\n'

		if type(pdb_input) == str:
			ncomp = 1
			if psf_input:
				input_str += f'coor @@{pdb_input}\n'
			else:
				pass

			# $file variable is still used by some CNS recipes, need refactoring!
			input_str += f'eval ($file=\"{pdb_input}\")\n'

		if type(pdb_input) == list or type(pdb_input) == tuple:
			ncomp = len(pdb_input)
			for pdb in pdb_input:
				input_str += f'coor @@{pdb}\n'

		chainsegs = PDB.identify_chainseg(pdb_input)

		ncomponents = len(chainsegs)

		input_str += f'eval ($ncomponents={ncomponents})\n'

		for i, segid in enumerate(chainsegs):
			input_str += f'eval ($prot_segid_mol{i+1}="{segid}")\n'

		try:
			ambig_fname = glob.glob('data/ambig.tbl')[0]
			input_str += f'eval ($ambig_fname="{ambig_fname}")\n'
		except IndexError:
			input_str += f'eval ($ambig_fname="")\n'

		try:
			unambig_fname = glob.glob('data/unambig.tbl')[0]
			input_str += f'eval ($unambig_fname="{unambig_fname}")\n'
		except IndexError:
			input_str += f'eval ($unambig_fname="")\n'

		try:
			hbond_fname = glob.glob('data/hbond.tbl')[0]
			input_str += f'eval ($hbond_fname="{hbond_fname}")\n'
		except IndexError:
			input_str += f'eval ($hbond_fname="")\n'

		try:
			dihe_fname = glob.glob('data/dihe.tbl')[0]
			input_str += f'eval ($dihe_fname="{dihe_fname}")\n'
		except IndexError:
			input_str += f'eval ($dihe_fname="")\n'

		try:
			tensor_fname = glob.glob('data/tensor.tbl')[0]
			input_str += f'eval ($tensor_tbl="{tensor_fname}")\n'
		except IndexError:
			input_str += f'eval ($dihe_fname="")\n'

		seed = random.randint(100, 999)
		input_str += f'eval ($seed={seed})\n'

		return input_str

	def prepare_output(self, pdb, psf, output_name):
		""" Tell the recipe wich should be the output file

		The output name must be UNIQUE
		"""
		output = '\n! Output structure\n'
		if output_name:
			if self.folder:
				if not os.path.isdir(self.folder):
					os.mkdir(self.folder)
			if psf:
				output += f"eval ($output_psf_filename= \"{self.folder}\" + \"{output_name}\" + \".psf\")\n"
			if pdb:
				output += f"eval ($output_pdb_filename= \"{self.folder}\" + \"{output_name}\" + \".pdb\")\n"
		else:
			print('+ ERROR: No output name defined for this recipe')
			exit()

		return output

	@staticmethod
	def load_protonation_state(prot_dic):
		""" Add protonation states to the recipe """
		protonation_header = '\n! Protonation states\n'

		for i, chain in enumerate(prot_dic):
			hise_l = [0] * 10
			hisd_l = [0] * 10
			hisd_counter = 0
			hise_counter = 0
			for res in prot_dic[chain]:
				state = prot_dic[chain][res].lower()
				if state == 'hise':
					hise_l[hise_counter] = res
					hise_counter += 1
				if state == 'hisd':
					hisd_l[hisd_counter] = res
					hisd_counter += 1

			hise_str = ''
			for e in [(i+1, c+1, r) for c, r in enumerate(hise_l)]:
				hise_str += f'eval ($toppar.hise_resid_{e[0]}_{e[1]} = {e[2]})\n'
			hisd_str = ''
			for e in [(i+1, c+1, r) for c, r in enumerate(hisd_l)]:
				hisd_str += f'eval ($toppar.hisd_resid_{e[0]}_{e[1]} = {e[2]})\n'

			protonation_header += hise_str
			protonation_header += hisd_str

		return protonation_header

	@staticmethod
	def load_ff_parameters(forcefield_parameters):
		""" Add force-field specific parameters to its appropriate places in the scoring recipe """
		ff_param_header = '\n! FF parameters\n'
		ff_param_header += 'parameter\n'
		ff_param_header += f'  @@{forcefield_parameters}\n'
		# ff_param_header += f'  @@{self.protein_param}\n'
		# ff_param_header += f'  @@{self.solvent_param}\n'
		# ff_param_header += f'  @@{self.ion_param}\n'
		# ff_param_header += f'  @@{self.ligand_param}\n'
		ff_param_header += 'end\n'

		return ff_param_header

	@staticmethod
	def load_ff_topology(forcefield_topology):
		""" Add force-field specific topology to its appropriate places in the scoring recipe """
		ff_top_header = '\n! Toplogy\n'
		ff_top_header += 'topology\n'
		ff_top_header += f'  @@{forcefield_topology}\n'
		# ff_top_header += f'  @@{self.protein_top}\n'
		# ff_top_header += f'  @@{self.solvent_top}\n'
		# ff_top_header += f'  @@{self.break_top}\n'
		# ff_top_header += f'  @@{self.ion_top}\n'
		# ff_top_header += f'  @@{self.ligand_top}\n'
		ff_top_header += 'end\n'

		return ff_top_header

	@staticmethod
	def load_link(mol_link):
		link_header = '\n! Link file\n'
		link_header += f'eval ($link_file = "{mol_link}" )\n'

		return link_header

	@staticmethod
	def load_trans_vectors(trans_vectors):
		trans_header = '\n! Translation vectors\n'
		for vector_id in trans_vectors:
			vector_file = trans_vectors[vector_id]
			trans_header += f'eval ($trans_vector_{vector_id} = "{vector_file}" )\n'

		return trans_header

	@staticmethod
	def load_tensor(tensor):
		tensor_header = '\n! Tensors'
		for tensor_id in tensor:
			tensor_file = tensor[tensor_id]
			tensor_header += f'eval (${tensor_id} = "{tensor_file}" )\n'

		return tensor_header

	@staticmethod
	def load_axis(axis):
		axis_header = '\n! Axis'
		for axis_id in axis:
			axis_file = axis[axis_id]
			axis_header = f'eval (${axis_id} = "{axis_file}" )\n'

		return axis_header

	@staticmethod
	def load_scatter(scatter_lib):
		scatter_header = '\n! Scatter lib'
		scatter_header += f'eval ($scatter_lib = "{scatter_lib}" )\n'

		return scatter_header
