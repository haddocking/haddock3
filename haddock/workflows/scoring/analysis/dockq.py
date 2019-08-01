import itertools
import subprocess
import haddock.workflows.scoring.config as config
from haddock.modules.structure.utils import PDB

segid2chain = PDB.segid2chain
identify_chains = PDB.identify_chains

dockq_exec = config.ini.get('third party', 'dockq_exe')


# TODO: Implement parallelism
def dockq(ref, pdb_f):
	print(f'+ {pdb_f}')
	irms = float('nan')
	lrms = float('nan')
	fnat = float('nan')
	capri = float('nan')
	dockq_score = float('nan')

	ref = segid2chain(ref)
	pdb_f = segid2chain(pdb_f)

	ref_comb = itertools.combinations(identify_chains(ref), 2)
	pdb_comb = itertools.combinations(identify_chains(pdb_f), 2)

	result_dic = {}
	for a, b in zip(pdb_comb, ref_comb):

		interface_name = '_'.join(a)

		cmd = f'{dockq_exec} {pdb_f} {ref} -model_chain1 {a[0]} {a[1]} -native_chain1 {b[0]} {b[1]}'

		p = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
		out = p.stdout.decode('utf-8').split('\n')

		for l in out:
			if 'Fnat' in l:
				fnat = float(l.split()[1])
			if 'iRMS' in l:
				irms = float(l.split()[1])
			if 'LRMS' in l:
				lrms = float(l.split()[1])
			if 'CAPRI ' in l:
				capri = l.split()[1]
			if 'DockQ ' in l and '*' not in l:
				dockq_score = float(l.split()[1])

		result_dic[f'{interface_name}_irms'] = irms
		result_dic[f'{interface_name}_lrms'] = lrms
		result_dic[f'{interface_name}_fnat'] = fnat
		result_dic[f'{interface_name}_capri'] = capri
		result_dic[f'{interface_name}_dockq'] = dockq_score

	return result_dic


#
# def define_interfaces(chain_str):
# 	interface_list = []
# 	if len(chain_str) >= 4:
# 		for inter_a in itertools.combinations(chain_str, 2):
# 			a = ''.join(inter_a)
# 			for inter_b in itertools.combinations(chain_str.replace(a, ''), 2):
# 				b = ''.join(inter_b)
# 				if not set(inter_a) & set(inter_b):
# 					if not (b, a) in interface_list:
# 						interface_list.append((a, b))
#
# 	elif len(chain_str) == 2:
# 		interface_list = list(itertools.combinations(chain_str, 2))
#
# 	elif len(chain_str) == 3:
# 		for inter_a in itertools.combinations(chain_str, 2):
# 			a = ''.join(inter_a)
# 			inter_b = set(chain_str) - set(inter_a)
# 			# inter_b = list(chain_str.replace(a, ''))
# 			b = ''.join(inter_b)
# 			if not set(inter_a) & set(inter_b):
# 				if not (b, a) in interface_list:
# 					interface_list.append((a, b))
#
# 	return interface_list
