import filecmp
import unittest
from haddock.modules.worker.distribution import JobCreator
from utils.files import get_full_path

data_path = get_full_path('test', 'test_data')


class TestJobCreator(unittest.TestCase):

	def setUp(self):
		self.jobs = JobCreator()

	def test_delegate(self):

		delegate_dic = self.jobs.delegate('recipe-as-string', ['structures/00000.pdb'])
		# local_path = os.path.dirname(__file__)
		dic = {0: ('jobs/00000.inp', 'out/00000.out')}

		k = list(delegate_dic.keys())[0]
		delegate_dic[k] = [('/'.join(a[0].split('/')[-2:]), '/'.join(a[1].split('/')[-2:])) for a in [delegate_dic[d] for d in delegate_dic]][0]

		self.assertEqual(delegate_dic, dic)
		self.assertTrue(filecmp.cmp(f'{data_path}/delegate-test.inp', delegate_dic[0][0]))


if __name__ == '__main__':
	unittest.main()
