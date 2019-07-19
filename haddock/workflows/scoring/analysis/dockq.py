import itertools
import subprocess

import haddock.workflows.scoring.config as config

from haddock.modules.structure.utils import PDB

segid2chain = PDB.segid2chain
identify_chains = PDB.identify_chains

# param_dic = load_parameters()
dockq_exec = config.ini.get('third party', 'dockq_exe')


# TODO: Implement parallelism
def dockq(ref, pdb_f):
	irms = float('nan')
	lrms = float('nan')
	fnat = float('nan')
	capri = float('nan')
	dockq_score = float('nan')

	segid2chain(ref)
	segid2chain(pdb_f)

	chain_l = ''.join(identify_chains(ref))
	interfaces = define_interfaces(chain_l)
	result_dic = {}
	for inter in interfaces:
		cmd = f'{dockq_exec} {ref} {pdb_f} -native_chain1 {inter[0][0]} {inter[0][1]} -model_chain1 {inter[0][0]} {inter[0][1]} ' \
			f'-native_chain2 {inter[1][0]} {inter[1][1]} -model_chain2 {inter[1][0]} {inter[1][1]}'

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
			if 'DockQ ' in l and not '*' in l:
				dockq_score = float(l.split()[1])

		result_dic['_'.join(inter)] = {'irms': irms, 'lrms': lrms, 'fnat': fnat, 'capri': capri, 'dockq_score': dockq_score}

	return result_dic


def define_interfaces(chain_str):
	interface_list = []
	if len(chain_str) >= 4:
		for inter_a in itertools.combinations(chain_str, 2):
			a = ''.join(inter_a)
			for inter_b in itertools.combinations(chain_str.replace(a, ''), 2):
				b = ''.join(inter_b)
				if not set(inter_a) & set(inter_b):
					if not (b, a) in interface_list:
						interface_list.append((a, b))

	elif len(chain_str) == 2:
		interface_list = list(itertools.combinations(chain_str, 2))

	elif len(chain_str) == 3:
		for inter_a in itertools.combinations(chain_str, 2):
			a = ''.join(inter_a)
			inter_b = set(chain_str) - set(inter_a)
			# inter_b = list(chain_str.replace(a, ''))
			b = ''.join(inter_b)
			if not set(inter_a) & set(inter_b):
				if not (b, a) in interface_list:
					interface_list.append((a, b))

	return interface_list
