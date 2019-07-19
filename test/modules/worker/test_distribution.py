import filecmp
import os
import unittest
from haddock.modules.worker.distribution import JobCreator

data_path = os.path.join(os.path.dirname(__file__), '../../test_data')


class TestStringMethods(unittest.TestCase):

	def setUp(self):
		self.jobs = JobCreator()

	def test_delegate(self):

		delegate_dic = self.jobs.delegate('recipe-as-string', ['structures/00000.pdb'])

		dic = {0: ('/home/rodrigo/haddock3/jobs/00000.inp', '/home/rodrigo/haddock3/out/00000.out')}

		self.assertEqual(delegate_dic, dic)
		self.assertTrue(filecmp.cmp(f'{data_path}/delegate-test.inp', delegate_dic[0][0]))


if __name__ == '__main__':
	unittest.main()