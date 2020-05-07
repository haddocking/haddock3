import itertools
import subprocess
from haddock.modules.structure.utils import PDB
from haddock.utils.files import get_full_path


def dockq(ref, pdb_f, dockq_exec):

	irms = .0
	lrms = .0
	fnat = .0
	dockq_score = .0
	capri = None
	order = None

	ref = PDB.fix_id(ref)
	pdb_f = PDB.fix_id(pdb_f)

	reference_chains = PDB.identify_chainseg(ref)
	pdb_chains = PDB.identify_chainseg(pdb_f)

	result_dic = {}

	if reference_chains != pdb_chains:
		print(f'+ WARNING: Skipping {pdb_f}, number of chains do not match. Expected {len(reference_chains)} found {len(pdb_chains)}')
		interface_name = ''
		result_dic[f'{interface_name}_irms'] = .0
		result_dic[f'{interface_name}_lrms'] = .0
		result_dic[f'{interface_name}_fnat'] = .0
		result_dic[f'{interface_name}_capri'] = ''
		result_dic[f'{interface_name}_dockq'] = .0
		result_dic[f'{interface_name}_order'] = ''

	else:

		for comb in itertools.combinations(pdb_chains, 2):

			# interface_name = ''.join(comb) + '-' + [e for e in pdb_chains if e not in comb][0]
			interface_name = ''.join(comb)

			cmd = f'{dockq_exec} {pdb_f} {ref} -native_chain1 {comb[0]} {comb[1]} -perm1'

			p = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
			out = p.stdout.decode('utf-8').split('\n')

			for l in out:
				if l.startswith('Best score'):
					output_l = l.split()
					order = f'{output_l[-4]}->{output_l[-1]}'
				elif l.startswith('Fnat'):
					fnat = float(l.split()[1])
				elif l.startswith('iRMS'):
					irms = float(l.split()[1])
				elif l.startswith('LRMS'):
					lrms = float(l.split()[1])
				elif l.startswith('DockQ_CAPRI'):
					capri = l.split()[1]
				elif l.startswith('DockQ'):
					dockq_score = float(l.split()[1])
				else:
					pass

			result_dic[f'{interface_name}_irms'] = irms
			result_dic[f'{interface_name}_lrms'] = lrms
			result_dic[f'{interface_name}_fnat'] = fnat
			result_dic[f'{interface_name}_capri'] = capri
			result_dic[f'{interface_name}_dockq'] = dockq_score
			result_dic[f'{interface_name}_order'] = order

	return result_dic

