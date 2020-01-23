# create the izone based on structural alignment
import multiprocessing
import sys
from _operator import itemgetter
import argparse
import shutil
import subprocess
import glob
import os
import random
import uuid
from itertools import groupby
# import operator
import tqdm


def split_chain(pdbf):
	# """ inpsired """ by https://github.com/JoaoRodrigues/pdb-tools (:

	prev_chain, chain_ids, chain_atoms = None, [], {}
	cur_chain = None
	for line in open(pdbf):
		if 'ATOM' in line[:4]:
			if prev_chain != line[21]:
				if not line[21] in chain_atoms:
					cur_chain = chain_atoms[line[21]] = []
				else:
					cur_chain = chain_atoms[line[21]]
				cur_chain.append(line)
				prev_chain = line[21]
				chain_ids.append(line[21])
			else:
				cur_chain.append(line)

	# Output chains to files
	pdb_dic = {}
	for c_id in chain_ids:
		name = str(uuid.uuid4())
		# name = pdbf.split('.pdb')[0] + '_' + c_id + '.pdb'
		pdb_dic[c_id] = name
		out = open(name, 'w')
		out.write(''.join(chain_atoms[c_id]))
		out.write('END\n')
		out.close()

	return pdb_dic


def load_seq(prot):
	aa_dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
			  'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
			  'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
			  'TYR': 'Y', 'VAL': 'V', 'DC': 'C',
			  'DA': 'A', 'DG': 'G', 'DT': 'T',
			  'ADE': 'A', 'THY':'T', 'GUA': 'G', 'CYT': 'C'}
	seq_dic = {}
	for l in open(prot):
		if 'ATOM' in l[:4]:
			chain = l[21]
			resnum = int(l[22:26])
			resname = l[17:20].split()[0]
			try:
				_ = seq_dic[chain]
			except KeyError:
				seq_dic[chain] = {}
			try:
				name = aa_dic[resname]
			except KeyError:
				name = 'X'
			seq_dic[chain][resnum] = name
	return seq_dic


def align(prota, protb, lovoalignexe='lovoalign'):

	pad = split_chain(prota)
	pbd = split_chain(protb)

	# check if chain ids match
	if pad.keys() != pbd.keys():
		print(f'+ ERROR: ChainIDs do not match! {prota} {protb}')
		exit()

	# make sure of numbering
	pa_seqdic = load_seq(prota)
	pb_seqdic = load_seq(protb)

	numbering_dic = {}
	for chain in pad.keys():
		numbering_dic[chain] = {}
		# cmd = '{} -p1 {} -p2 {} -c1 {} -c2 {}'.format(lovoalignexe, pad[chain], pbd[chain], chain, chain)
		cmd = '{} -p1 {} -p2 {} -c1 {} -c2 {} -seqnum'.format(lovoalignexe, pad[chain], pbd[chain], chain, chain)

		out = subprocess.getoutput(cmd).split('\n')

		aln_start = [out.index(e) for e in out if 'SEQUENCE ALIGNMENT' in e][0] + 2
		aln_end = [i - 2 for i, k in enumerate(out[aln_start:]) if 'FINAL' in k][0] + aln_start
		aln_l = out[aln_start:aln_end]

		aln = [aln_l[i:i + 3][:2] for i in range(0, len(aln_l), 3)]

		idx_a = 0
		idx_b = 0

		for e in aln:
			a, b = e

			a = a.split()
			b = b.split()

			a_seq = a[1]
			b_seq = b[1]

			reslist_a = list(pa_seqdic[chain].keys())
			reslist_b = list(pb_seqdic[chain].keys())

			res_a = None
			res_b = None

			for a_aa, b_aa in zip(a_seq, b_seq):
				if a_aa != '-':
					res_a = reslist_a[idx_a]
					idx_a += 1

				if b_aa != '-':
					res_b = reslist_b[idx_b]
					idx_b += 1

				if a_aa != '-' and b_aa != '-':
					numbering_dic[chain][res_a] = res_b
	to_be_removed = ' '.join(list(pad.values()) + list(pbd.values()))
	os.system(f'rm {to_be_removed}')

	return numbering_dic


def output_renumbered(prot, numbering_dic):
	#
	renumb_pdb_l = []
	for l in open(prot):
		if 'ATOM' in l[:4]:
			chain = l[21]
			resnum = int(l[22:26])
			new_res = numbering_dic[chain][resnum]
			n_l = l[:22] + '{:>4}'.format(new_res) + ' ' + l[26:]
			renumb_pdb_l.append(n_l)

	outputf = prot.replace('.pdb', '-renum.pdb')
	out = open(outputf, 'w')
	out.write(''.join(renumb_pdb_l))
	out.close()
	print(prot)


def run_contacts(pdbf, cutoff):
	cmd = 'contact {} {}'.format(pdbf, cutoff)
	out = subprocess.getoutput(cmd).split('\n')
	return out


def identify_inteface(pdbf, cutoff):
	contacts_l = run_contacts(pdbf, cutoff)

	interface_dic = {}
	for l in contacts_l:
		resnum_a, chain_a, atom_a, resnum_b, chain_b, atom_b, distance = l.split()
		resnum_a = int(resnum_a)
		resnum_b = int(resnum_b)

		# One way
		try:
			_ = interface_dic[chain_a]
		except KeyError:
			interface_dic[chain_a] = {}
		try:
			_ = interface_dic[chain_a][chain_b]
		except KeyError:
			interface_dic[chain_a][chain_b] = []

		if float(distance) <= cutoff:
			if resnum_a not in interface_dic[chain_a][chain_b]:
				interface_dic[chain_a][chain_b].append(resnum_a)

		# other way
		try:
			_ = interface_dic[chain_b]
		except KeyError:
			interface_dic[chain_b] = {}
		try:
			_ = interface_dic[chain_b][chain_a]
		except KeyError:
			interface_dic[chain_b][chain_a] = []

		if float(distance) <= cutoff:
			if resnum_b not in interface_dic[chain_b][chain_a]:
				interface_dic[chain_b][chain_a].append(resnum_b)

	# sort residue lists
	ninterface_dic = {}
	for a in interface_dic:
		ninterface_dic[a] = {}
		for b in interface_dic[a]:
			reslist = interface_dic[a][b]
			reslist.sort()
			ninterface_dic[a][b] = reslist

	return ninterface_dic


def get_range(data):
	ranges = []
	for k, g in groupby(enumerate(data), lambda x: x[0] - x[1]):
		group = (map(itemgetter(1), g))
		group = list(map(int, group))
		ranges.append((group[0], group[-1]))
	return ranges


def retrieve_izone(c_dic, numbering_dic):
	# based on the reference interface, create izone
	izone_l = []
	for chain in c_dic:
		ref_dic = {}
		for bound_res in list(c_dic[chain].items())[0][1]:
			try:
				ub = numbering_dic[chain][bound_res]
				ref_dic[bound_res] = ub
			except KeyError:
				pass

		for bound_range in get_range(ref_dic.keys()):
			unbound_res_l = []
			for bound_res in range(bound_range[0], bound_range[1] + 1):
				unbound_res_l.append(ref_dic[bound_res])

			for unbound_range in get_range(unbound_res_l):
				bound_res_l = []
				for unbound_res in range(unbound_range[0], unbound_range[1] + 1):
					bound_res_l.append(list(ref_dic.keys())[list(ref_dic.values()).index(unbound_res)])

				range_a = get_range(bound_res_l)[0]  # bound
				range_b = unbound_range

				izone_str = 'ZONE %s%i-%s%i:%s%i-%s%i' % (
					chain, range_a[0], chain, range_a[1], chain, range_b[0], chain, range_b[1])
				izone_l.append(izone_str)

	return izone_l


def run_profit(cmd):
	return subprocess.getoutput('echo "{}" | profit'.format(cmd)).split('\n')


def calc_irmsd(prot_a, prot_b, atoms, numbering_dic):

	contact_dic_a = identify_inteface(prot_a, 10.0)
	izone_l = retrieve_izone(contact_dic_a, numbering_dic)

	# cmd = 'refe {}\nmobi {}\n{}\natoms {}\nfit\nzone clear'.format(prot_a, prot_b, '\n'.join(izone_l), atoms)
	cmd = 'refe %s\nmobi %s\nATOMS %s\nZONE CLEAR\n%s\nstatus\nFIT\nquit' % (prot_a, prot_b, atoms, '\n'.join(izone_l))
	# open('irmsd.dbg', 'w').write(cmd)
	# open('izone', 'w').write('\n'.join(izone_l))

	out = run_profit(cmd)

	irmsd = None
	try:
		irmsd = float(out[-2].split()[-1])
	except KeyError:
		print('Something went wrong when running PROFIT, check irmsd.dbg')
		exit()

	return irmsd


def calc_fnat(pa, pb, numbering_dic, cutoff=5.0):
	con_a = run_contacts(pa, cutoff)
	con_b = run_contacts(pb, cutoff)

	a_con_l = []
	b_con_l = []
	for e in con_a:
		resnum_x, chain_x, _, resnum_y, chain_y, _, _ = e.split()
		try:
			resnum_x = str(numbering_dic[chain_x][int(resnum_x)])
			resnum_y = str(numbering_dic[chain_y][int(resnum_y)])
		except KeyError:
			# one of the residues present in this contact was not matched to the target
			continue

		a_con_l.append((resnum_x, resnum_y))
	a_con_l = set(a_con_l)

	for e in con_b:
		resnum_x, _, _, resnum_y, _, _, _ = e.split()

		b_con_l.append((resnum_x, resnum_y))
	b_con_l = set(b_con_l)


	try:
		fnat = len(a_con_l & b_con_l) / len(a_con_l)
	except ZeroDivisionError:
		# No contacts were matched
		fnat = .0

	return fnat


def calc_lrmsd(prot_a, prot_b, numbering_dic, atoms):
	lrmsd = None

	# receptor = first chain
	# retrieve_lzone(numbering_dic)
	chain_l = list(numbering_dic.keys())
	chain_l.sort()
	receptor_chain = chain_l[0]
	ligand_zone = {}
	for chain in numbering_dic:
		ligand_zone[chain] = []
		for bound_range in get_range(numbering_dic[chain]):
			unbound_res_l = []
			for bound_res in range(bound_range[0], bound_range[1] + 1):
				unbound_res_l.append(numbering_dic[chain][bound_res])

			for unbound_range in get_range(unbound_res_l):
				bound_res_l = []
				for unbound_res in range(unbound_range[0], unbound_range[1] + 1):
					bound_res_l.append(list(numbering_dic[chain].keys())[list(numbering_dic[chain].values()).index(unbound_res)])

				range_a = get_range(bound_res_l)[0]  # bound
				range_b = unbound_range

				lzone_str = 'ZONE %s%i-%s%i:%s%i-%s%i' % (
					chain, range_a[0], chain, range_a[1], chain, range_b[0], chain, range_b[1])
				ligand_zone[chain].append(lzone_str)
				# lzone_l.append(lzone_str)


		# for ref_res in numbering_dic[chain]:
		#     target_res = numbering_dic[chain][ref_res]
		#     lzone = 'ZONE %s%s-%s%i:%s%i-%s%i' % (chain, ref_res, chain, ref_res, chain, target_res, chain, target_res)
		#     ligand_zone[chain].append(lzone)

	lzone_str_out = '\n'.join(ligand_zone[receptor_chain])
	lzone_str_out += '\n'

	cmd = 'refe {}'.format(prot_a)
	cmd += '\n'
	cmd += 'mobi {}'.format(prot_b)
	cmd += '\n'
	cmd += '\n'.join(ligand_zone[receptor_chain])
	cmd += '\n'
	cmd += 'atoms {}'.format(atoms)
	cmd += '\n'
	cmd += 'fit'
	cmd += '\n'
	for ligand in ligand_zone:
		if ligand != receptor_chain:
			l_tbw = ''
			for zone in ligand_zone[ligand]:
				l_tbw += ' R%s\n' % zone
				lzone_str_out += ' R%s\n' % zone
			cmd += l_tbw[1:]
			cmd += '\n'
	cmd += 'ZONE CLEAR'
	cmd += '\n'
	cmd += 'quit'

	# open('lrmsd.dbg', 'w').write(cmd)
	# open('lzone','w').write(lzone_str_out)

	out = run_profit(cmd)
	try:
		lrmsd = float([e for e in out if 'RMS' in e][-1].split()[-1])
	except KeyError:
		print('Something went wrong when running PROFIT, check irmsd.dbg')
		exit()

	return lrmsd


def clean(prot_a, prot_b):
	pdb_l = [glob.glob('{}_*'.format(e.split('.pdb')[0])) for e in [prot_a, prot_b]]
	pdb_l = [x for xs in pdb_l for x in xs] + glob.glob('*flatnum*')
	for p in set(pdb_l):
		os.system(f'rm {p}')


def scramble_prot(prot):
	# scramble protein numbering for debug purposes
	seq = load_seq(prot)
	chain_resdic = dict([(c, list(seq[c].keys())) for c in seq])
	new_resdic = {}
	for chain in chain_resdic:
		new_resdic[chain] = {}
		reslist = chain_resdic[chain]
		shuffled_reslist = reslist.copy()
		random.shuffle(shuffled_reslist)
		new_resdic[chain] = dict(zip(reslist, shuffled_reslist))

	scrambled_prot = []
	for l in open(prot):
		if 'ATOM' in l[:4]:
			chain = l[21]
			resnum = int(l[22:26])
			new_resnum = new_resdic[chain][resnum]
			n_l = l[:22] + '{0:>4}'.format(new_resnum) + l[26:]
			scrambled_prot.append(n_l)

	outname = '{}_scramb.pdb'.format(prot.split('.pdb')[0])
	out = open(outname, 'w')
	out.write(''.join(scrambled_prot))
	out.close()

	return outname


def flatten_numbers(prot):
	# look for residue insertions
	# 10, <10A, 10B, 10C>, 12

	# create a numbering dictionary taking into account insertions
	resdic = {}
	chain_l = []
	incr = None
	for l in open(prot):
		if 'ATOM' in l[:4] and 'CA' in l[12:16]:

			chain = l[21]

			if not chain in chain_l:
				incr = 0
				chain_l.append(chain)

			try:
				_ = resdic[chain]
			except KeyError:
				resdic[chain] = {}

			ori_resnum = int(l[22:26])
			icode = l[26]

			if not icode.isspace():
				incr += 1

			resnum = ori_resnum + incr
			resdic[chain][(ori_resnum, icode)] = resnum

	flatf_l = []
	for l in open(prot):
		if 'ATOM' in l[:4]:
			chain = l[21]
			resnum = int(l[22:26])
			icode = l[26]
			new_res = resdic[chain][(resnum, icode)]
			n_l = l[:22] + '{:>4}'.format(new_res) + ' ' + l[27:]
			flatf_l.append(n_l)

	outputf = prot.split('.pdb')[0] + '_flatnum.pdb'
	out = open(outputf,'w')
	out.write(''.join(flatf_l))
	out.close()

	return outputf


def get_pdb_list(file_list):
	# when file.list changes, edit this function!
	# currently handling: "PREVIT:complex_it0_000009.pdb"  { -7.739781000000001 }
	path = '/'.join(file_list.split('/')[:-1])
	result = []
	with open(file_list) as f:
		for l in f.readlines():
			pdb = l.split()[0].split(':')[1].split('"')[0]
			if not os.path.isfile(pdb):
				if os.path.isfile(f'{path}/{pdb}'):
					pdb = f'{path}/{pdb}'
				else:
					print('PDB not found, check under the hood')
					exit()
			haddock_score = float(l.split()[-2])
			result.append((pdb, haddock_score))
	return result


def main():
	error_check = False
	# lovoalign_exe = '/home/rodrigo/software/lovoalign'
	for exe in ['profit', 'contact', 'lovoalign']:
		if not shutil.which(exe):
			print('ERROR: {} not found in $PATH'.format(exe))
			error_check = True

	if error_check:
		exit()

	parser = argparse.ArgumentParser()
	parser.add_argument("prot_a", help="Reference protein")
	# parser.add_argument("prot_b", help="Target protein")
	parser.add_argument("pdb_file_list", help="file.list")
	parser.add_argument("output_file", help="output filename")

	# parser.add_argument("--cg", help="Use CG beads", action="store_true", default=False)
	# parser.add_argument("--flat", help="Flatten numbering, use this if your protein has insertions ex: 10, 10A, 10B, 11")
	parser.add_argument("--nproc", help="How many processors should be used", type=int)

	args = parser.parse_args()

	if args.nproc:
		nproc = args.nproc
		print(f'+ Using {nproc} processors')
	else:
		nproc = multiprocessing.cpu_count()
		print(f'+ Using {nproc} processors')

	# if args.cg:
	# 	atoms = 'BB*'
	# else:
	# 	atoms = 'CA,N,C,O'

	pa = args.prot_a

	result_list = get_pdb_list(args.pdb_file_list)

	arg_list = []
	for element in result_list:
		pb, haddock_score = element
		arg_list.append((pa, pb, haddock_score))

	# result = []
	# with multiprocessing.Pool(processes=nproc) as pool:
	# 	for i, out in enumerate(pool.imap(execute, arg_list)):
	# 		per = (i / len(arg_list)) * 100
	# 		sys.stderr.write(f'\r++ Progress: {per:.2f} %')
	# 		result.append(out)

	with multiprocessing.Pool(processes=nproc) as pool:
		result = list(tqdm.tqdm(pool.imap(execute, arg_list), total=len(arg_list)))

	with open(args.output_file, 'w') as f:
		f.write('pdb\tirmsd\tlrmsd\tfnat\n')
		for l in result:
			f.write(l + '\n')
	f.close()


def execute(args):
	pa, pb, haddock_score = args
	num_dic = align(pa, pb)

	irmsd = calc_irmsd(pa, pb, 'CA,N,C,O', num_dic)
	fnat = calc_fnat(pa, pb, num_dic)
	lrmsd = calc_lrmsd(pa, pb, num_dic, 'CA,N,C,O')

	return f'{pb}\t{haddock_score:.3f}\t{irmsd:.4f}\t{lrmsd:.4f}\t{fnat:.4f}'


if __name__ == '__main__':
	main()
