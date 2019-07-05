from haddock.workflows.scoring.config import load_parameters


class RecipeGenerator:

	def __init__(self):
		self.input_header = None
		self.ff_param_header = None
		self.ff_top_header = None
		self.scoring_header = None
		self.link_header = None
		self.header = None
		self.param_dic = {}
		self.recipe = None

	def initialize(self):
		# with open(param_file) as handle:
		# 	self.param_dic = json.loads(handle.read())
		self.param_dic = load_parameters()
		self.generate_input_header()
		self.load_ff_parameters()
		self.load_ff_topology()
		self.load_scoring_parameters()
		self.load_link()
		self.combine_header()

	def generate_input_header(self):
		""" Create input header with number of components and its segids """
		ncomponents = self.param_dic['components']['number']
		self.input_header = f'evaluate($data.ncomponents={ncomponents})\n'

		for i, seg in enumerate(self.param_dic['components']['segids']):
			self.input_header += f'evaluate($Toppar.prot_segid_{i} = "{seg}")\n'

	def load_ff_parameters(self):
		""" Add force-field specific parameters to its appropriate places in the scoring recipe """
		self.ff_param_header = 'parameter\n'
		for param in self.param_dic['ff-parameters']:
			v = self.param_dic['ff-parameters'][param]
			self.ff_param_header += f'  @@{v}\n'
		self.ff_param_header += 'end\n'

	def load_ff_topology(self):
		""" Add force-field specific topology to its appropriate places in the scoring recipe """
		self.ff_top_header = 'topology\n'
		for top in self.param_dic['topology']:
			v = self.param_dic['topology'][top]
			self.ff_top_header += f'  @@{v}\n'
		self.ff_top_header += 'end\n'

	def load_scoring_parameters(self):
		self.scoring_header = ''
		for flag in self.param_dic['scoring-parameters']['flags']:
			v = str(self.param_dic['scoring-parameters']['flags'][flag]).upper()
			self.scoring_header += f'evaluate($Data.flags.{flag} = {v})\n'

		for value in self.param_dic['scoring-parameters']['values']:
			v = self.param_dic['scoring-parameters']['values'][value]
			self.scoring_header += f'evaluate(${value}={v})\n'

	def load_link(self):
		v = self.param_dic['link']
		self.link_header = f'evaluate ($link_file = "{v}" )'

	def combine_header(self):
		""" Create the final header to be used by the scoring recipe """
		self.header = self.input_header + self.ff_param_header + self.ff_top_header + self.scoring_header + self.link_header

	def load(self):
		""" Load a CNS script and finalize the recipe """
		recipe_f = self.param_dic['input']['recipe']
		with open(recipe_f) as f:
			self.recipe = self.header + ''.join(f.readlines())
		f.close()

