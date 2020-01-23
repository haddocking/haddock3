import unittest
from haddock.modules.worker.recipe import RecipeComposer
from haddock.utils.files import get_full_path

data_path = get_full_path('test', 'data')


class TestRecipeComposer(unittest.TestCase):
	def setUp(self):

		default = {"params": {"log_level": "verbose"}}
		custom = {"params": {"log_level": "quiet"}}

		self.Recipe = RecipeComposer(recipe_file='rigid-body.cns', default_params=default, custom_params=custom)

	def test_identify_modules(self):

		recipe_f = f'{data_path}/test_recipe.inp'
		module_list = self.Recipe.identify_modules(recipe_f)
		expected_module_list = ['auto-his.cns']

		self.assertEqual(module_list, expected_module_list)

	def test_list_dependencies(self):
		recipe_f = f'{data_path}/test_recipe.inp'
		dependency_list = self.Recipe.list_dependencies(recipe_f)
		expected_dependency_list = ['auto-his.cns']

		self.assertEqual(dependency_list, expected_dependency_list)

	def test_compose(self):

		recipe_str = self.Recipe.compose()

		# ENHANCEMENT: Test this in a more complete manner
		segment_10_20 = 'ers\neval ('
		segment_30_40 = '="quiet")\n'
		segment_60_70 = 'e core doc'
		segment_90_100 = 'CK perform'

		self.assertEqual(recipe_str[10:20], segment_10_20)
		self.assertEqual(recipe_str[30:40], segment_30_40)
		self.assertEqual(recipe_str[60:70], segment_60_70)
		self.assertEqual(recipe_str[90:100], segment_90_100)

	def test_load_recipe_params(self):

		param = self.Recipe.load_recipe_params()
		expected_param = '\n! Parameters\neval ($log_level="quiet")\n'

		self.assertEqual(param, expected_param)


if __name__ == '__main__':
	unittest.main()
