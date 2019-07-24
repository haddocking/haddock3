""" Test the Input generator """
import unittest
from unittest.mock import patch
from utils.files import get_full_path

data_path = get_full_path('test', 'test_data')
param_path = get_full_path('haddock', 'toppar')

with patch("sys.argv", ['', f'{data_path}/scoring.json']):
	from haddock.modules.cns.input import HeaderComposer, RecipeComposer, RecipeGenerator


class TestHeaderComposer(unittest.TestCase):

	def setUp(self):
		self.hc = HeaderComposer()

	def test_load_ff_parameters(self):
		param = self.hc.load_ff_parameters()
		param_header_str = f'parameter\n  @@{param_path}/protein-allhdg5-4.param\n  @@{param_path}/water-allhdg5-4.param\n  @@{param_path}/ion.param\n  @@{param_path}/ligand.param\nend\n'

		self.assertEqual(param, param_header_str)

	def test_load_ff_topology(self):
		top = self.hc.load_ff_topology()

		top_header_str = f'topology\n  @@{param_path}/protein-allhdg5-4.top\n  @@{param_path}/water-allhdg5-4.top\n  @@{param_path}/protein_break.top\n  @@{param_path}/ion.top\n  @@{param_path}/ligand.top\nend\n'

		self.assertEqual(top, top_header_str)

	def test_load_scoring_parameters(self):
		scoring = self.hc.load_scoring_parameters()

		scoring_header_str = 'evaluate($Data.flags.dihed = FALSE)\nevaluate($Data.flags.sani = FALSE)\nevaluate($Data.flags.coup = FALSE)\nevaluate($Data.flags.vean = FALSE)\nevaluate($Data.flags.cdih = FALSE)\nevaluate($Data.flags.noe = TRUE)\nevaluate($Data.flags.sym = FALSE)\nevaluate($Data.flags.ncs = FALSE)\nevaluate($Data.flags.noecv = FALSE)\nevaluate($Data.flags.auto_break = TRUE)\nevaluate($Data.flags.waterdock = TRUE)\nevaluate($break_cutoff=2.5)\nevaluate($hydrogen_build=all)\nevaluate($disulphide_dist=3)\nevaluate($log_level=verbose)\nevaluate($epsilon=1)\n'

		self.assertEqual(scoring, scoring_header_str)

	def test_load_link(self):
		link = self.hc.load_link()
		link_header_str = f'evaluate ($link_file = "{param_path}/protein-allhdg5-4-noter.link" )'

		self.assertEqual(link, link_header_str)

	def test_protonation_dic(self):
		dummy_dic = {'A': {42: 'HISE'}, 'B': {1: 'HISD', 2: 'HISE'}}
		protonation_str = self.hc.load_protonation_state(dummy_dic)

		expected_str = 'evaluate ($toppar.hise_resid_1_1 = 42)\nevaluate ($toppar.hisd_resid_2_1 = 1)\nevaluate ($toppar.hise_resid_2_2 = 2)\n'
		self.assertEqual(protonation_str, expected_str)

	def test_create_header(self):
		dummy_dic = {}
		header = self.hc.create_header(dummy_dic)

		self.assertEqual(header.split()[0], 'parameter')
		self.assertEqual(header.split()[6], 'topology')
		self.assertEqual(header.split()[-4], '($link_file')


class TestRecipeComposer(unittest.TestCase):

	def setUp(self):
		self.rc = RecipeComposer()

	def test_compose(self):
		body_str = self.rc.compose()

		with open(f'{data_path}/composed-recipe.cns') as f:
			composed = f.readlines()
		f.close()

		composed = ''.join(composed)

		self.assertEqual(body_str, composed)

	def test_identify_modules(self):
		# is it possible to test static methods?
		module_list = self.rc.identify_modules(f'{data_path}/dependency.cns')

		self.assertEqual(module_list, ['def_solv_param.cns'])

	def test_list_dependencies(self):
		dependency_list = self.rc.list_dependencies(f'{data_path}/dependency.cns')

		self.assertEqual(dependency_list, ['def_solv_param.cns'])


class TestRecipeGenerator(unittest.TestCase):

	def setUp(self):
		self.rg = RecipeGenerator()

	def test_generate(self):

		dummy_dic = {}
		generated_recipe_str = self.rg.generate(dummy_dic)

		generated_recipe_list = generated_recipe_str.split('\n')

		self.assertEqual(generated_recipe_list[0], 'parameter')
		self.assertEqual(generated_recipe_list[6], 'topology')
		self.assertEqual(generated_recipe_list[28], 'evaluate($epsilon=1)')


if __name__ == '__main__':
	unittest.main()

