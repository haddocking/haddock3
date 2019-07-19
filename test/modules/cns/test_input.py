""" Test the Input generator """
import unittest
from unittest.mock import patch
from utils.files import get_full_path

data_path = get_full_path('test', 'test_data')
param_path = get_full_path('haddock', 'toppar')

with patch("sys.argv", ['', f'{data_path}/scoring.json']):
	from haddock.modules.cns.input import HeaderComposer, RecipeComposer, RecipeGenerator


class TestStringMethods(unittest.TestCase):

	def setUp(self):
		self.hc = HeaderComposer()
		self.rc = RecipeComposer()
		self.rg = RecipeGenerator()

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

	def test_create_header(self):

		header = self.hc.create_header()

		self.assertEqual(header.split()[0], 'parameter')
		self.assertEqual(header.split()[6], 'topology')
		self.assertEqual(header.split()[-4], '($link_file')

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

	def test_generate(self):
		# FIXME: Refactor this test
		pass


if __name__ == '__main__':
	unittest.main()

