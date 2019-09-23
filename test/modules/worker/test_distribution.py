import unittest
import os
from haddock.modules.worker.distribution import JobCreator


class TestJobCreator(unittest.TestCase):
	def setUp(self):
		self.Job = JobCreator(job_id='test', job_folder='test_f')

	def test_delegate(self):

		delegated_dic = self.Job.delegate(job_num=1, input_file_str='')
		expected_dic = {1: ('test_f/test_000001.inp', 'test_f/test_000001.out')}

		self.assertEqual(delegated_dic, expected_dic)
		self.assertTrue(os.path.isfile('test_f/test_000001.inp'))

		os.remove('test_f/test_000001.inp')
		os.rmdir('test_f')


if __name__ == '__main__':
	unittest.main()
