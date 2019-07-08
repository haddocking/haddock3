import glob
# import operator
# import itertools
import re
import os

from haddock.workflows.scoring.analysis.contact import Contacts
from haddock.workflows.scoring.analysis.dfire import dfire
from haddock.workflows.scoring.analysis.dockq import dockq
from haddock.workflows.scoring.analysis.fastcontact import fastcontact
# from haddock.workflows.scoring.src import calc_fcc_matrix, cluster_fcc
# from haddock.workflows.scoring.analysis.profit import profit_rmsd
# from haddock.workflows.scoring.src import calc_fcc_matrix
from haddock.workflows.scoring.config import load_parameters


class Ana:

	def __init__(self):
		self.structure_dic = {}
		self.contact_filelist = []
		self.structure_haddockscore_list = []
		self.fcc_matrix_f = 'fcc.matrix'
		self.fcc_matrix = []
		self.con_list = None
		self.param_dic = load_parameters()

	def retrieve_structures(self):
		""" Retrieve structures that have been trough the CNS recipe """

		conv_l = glob.glob('structures/*conv.pdb')
		for pdb in conv_l:
			self.structure_dic[pdb] = {}

	def extract_energies(self):
		""" Extract energies from the header of the PDB file, according to HADDOCK formatting """

		total = float('nan')
		bonds = float('nan')
		angles = float('nan')
		improper = float('nan')
		dihe = float('nan')
		vdw = float('nan')
		elec = float('nan')
		air = float('nan')
		cdih = float('nan')
		coup = float('nan')
		sani = float('nan')
		vean = float('nan')
		dani = float('nan')
		desolv = float('nan')
		bsa = float('nan')

		for pdb in self.structure_dic:

			vdw_elec_air_regex = r"\s(\-?\d*\.?\d{1,}(E-\d{1,}|)|0\b|\$\w{1,})"
			desolv_regex = r"(\-?\d*\.?\d*)$"
			bsa_regex = r"(\-?\d*\.?\d*)$"

			f = open(pdb, 'r')
			for line in f:

				if 'REMARK energies' in line:
					# dirty fix to account for 8.754077E-02 notation and the eventual $DANI or $NOE
					energy_v = re.findall(vdw_elec_air_regex, line)

					temp_v = []
					for v in energy_v:
						v = v[0]
						try:
							v = float(v)
						except ValueError:
							v = .0
						temp_v.append(v)
					energy_v = temp_v
					###
					total, bonds, angles, improper, dihe, vdw, elec, air, cdih, coup, sani, vean, dani = energy_v
					vdw = float(vdw)
					elec = float(elec)
					air = float(elec)

				if 'REMARK Desolvation' in line:
					# print(line)
					desolv = float(re.findall(desolv_regex, line)[0])

				if 'REMARK buried surface area' in line:
					# print(line)
					bsa = float(re.findall(bsa_regex, line)[0])
					break

			f.close()

			self.structure_dic[pdb]['total'] = total
			self.structure_dic[pdb]['bonds'] = bonds
			self.structure_dic[pdb]['angles'] = angles
			self.structure_dic[pdb]['improper'] = improper
			self.structure_dic[pdb]['dihe'] = dihe
			self.structure_dic[pdb]['vdw'] = vdw
			self.structure_dic[pdb]['elec'] = elec
			self.structure_dic[pdb]['air'] = air
			self.structure_dic[pdb]['cdih'] = cdih
			self.structure_dic[pdb]['coup'] = coup
			self.structure_dic[pdb]['sani'] = sani
			self.structure_dic[pdb]['vean'] = vean
			self.structure_dic[pdb]['dani'] = dani
			self.structure_dic[pdb]['desolv'] = desolv
			self.structure_dic[pdb]['bsa'] = bsa

	def calculate_haddock_score(self):
		""" Calculate the HADDOCK Score of the PDB file using its appropriate weight """

		print(f'+ Calculating HADDOCK score for {len(self.structure_dic)} structures')

		for pdb in self.structure_dic:
			vdw = self.structure_dic[pdb]['vdw']
			elec = self.structure_dic[pdb]['elec']
			desolv = self.structure_dic[pdb]['desolv']
			air = self.structure_dic[pdb]['air']
			bsa = self.structure_dic[pdb]['bsa']

			# TODO: check if these weights are correct
			haddock_score = 1.0 * vdw + 0.2 * elec + 1.0 * desolv + 0.1 * air - 0.0 * bsa

			self.structure_dic[pdb]['haddock-score'] = haddock_score

	def run_fastcontact(self):
		""" Run fastcontact on all scored PDBs """

		print('+ Running FASTCONTACT')

		for pdb in self.structure_dic:
			fast_elec, fast_desol = fastcontact(pdb)
			self.structure_dic[pdb]['fastelec'] = fast_elec
			self.structure_dic[pdb]['fastdesol'] = fast_desol

	def run_dfire(self):
		""" Run dfire on all scored PDBs """

		print('+ Running DFIRE')

		for pdb in self.structure_dic:
			d_binding, d_score = dfire(pdb)
			self.structure_dic[pdb]['dfire-ebinding'] = d_binding
			self.structure_dic[pdb]['dfire-score'] = d_score

	def cluster(self, cutoff, strictness=0.75, threshold=4):
		""" Cluster scored models using FCC, output a sorted text file containing clusters and mean scores """

		print(f'+ Clustering with cutoff: {cutoff} and threshold: {threshold}')

		self.calculate_contacts()
		self.calc_fcc_matrix()

		cutoff = float(cutoff)
		partner_cutoff = float(cutoff) * float(strictness)

		elements = {}

		for e in self.fcc_matrix:
			ref, mobi, dRM, dMR = e.split()
			ref = int(ref)
			mobi = int(mobi)
			dRM = float(dRM)
			dMR = float(dMR)

			# Create or Retrieve Elements
			if ref not in elements:
				r = Element(ref)
				elements[ref] = r
			else:
				r = elements[ref]

			if mobi not in elements:
				m = Element(mobi)
				elements[mobi] = m
			else:
				m = elements[mobi]

			# Assign neighbors
			if dRM >= cutoff and dMR >= partner_cutoff:
				r.add_neighbor(m)
			if dMR >= cutoff and dRM >= partner_cutoff:
				m.add_neighbor(r)

		clusters = []
		threshold -= 1  # Account for center
		# ep = element_pool
		ep = elements
		cn = 1  # Cluster Number
		while 1:
			# Clusterable elements
			ce = [e for e in ep if not ep[e].cluster]
			if not ce:  # No more elements to cluster
				break

			# Select Cluster Center
			# Element with largest neighbor list
			ctr_nlist, ctr = sorted([(len([se for se in ep[e].neighbors if not se.cluster]), e) for e in ce])[-1]

			# Cluster until length of remaining elements lists are above threshold
			if ctr_nlist < threshold:
				break

			# Create Cluster
			c = Cluster(cn, ep[ctr])
			cn += 1
			clusters.append(c)

		cluster_dic = {}
		for c in clusters:
			clustered_models_list = []
			for model in [m.name for m in c.members]:
				model_str = '0' * (6 - len(str(model - 1))) + str(model - 1)
				name = f'structures/{model_str}_conv.pdb'
				haddock_score = self.structure_dic[name]['haddock-score']
				clustered_models_list.append((model, haddock_score))

			# sort cluster elements by haddock score
			model_list = sorted(clustered_models_list, key=lambda x: x[1])
			score_list = [e[1] for e in clustered_models_list]
			mean_score = sum(score_list) / len(score_list)
			# top4_mean_score = sum(score_list[:4]) / len(score_list[:4])
			tbw = f"Cluster {c.name} -> ({c.center.name}) "
			for m in model_list:
				tbw += f'{m[0]} '
			tbw += '\n'
			cluster_dic[c.name] = (tbw, mean_score)

		# sort clusters by mean haddock score
		sorted_cluster_list = sorted([(x, cluster_dic[x][1]) for x in list(cluster_dic.keys())], key=lambda x: x[1])

		# output
		with open(f'cluster_{cutoff}_{threshold + 1}.out', 'w') as out:
			for c in sorted_cluster_list:
				cluster_name, cluster_mean = c
				tbw_l = cluster_dic[cluster_name][0].split('->')
				tbw = f'{tbw_l[0]}[{cluster_mean:.3f}] ->{tbw_l[1]}'
				out.write(tbw)
		out.close()

	def calculate_contacts(self):
		""" Calculate atomic contacts """

		pdb_list = list(self.structure_dic.keys())

		con = Contacts()
		con.calculate_contacts(pdb_list, cutoff=5.0)

		self.con_list = con.contact_file_list

	def calc_fcc_matrix(self):
		""" Calculate the FCC matrix (extracted and adapted from calc_fcc_matrix.py """

		if not os.path.isfile(self.fcc_matrix_f) and not self.fcc_matrix:

			# print('+ Creating FCC matrix')

			self.con_list.sort()  # very important!
			contacts = [set([int(l) for l in open(f)]) for f in self.con_list if f.strip()]

			# get the pairwise matrix
			for i in range(len(contacts)):
				for j in range(i+1, len(contacts)):
					con_a = contacts[i]
					con_b = contacts[j]

					cc = float(len(con_a.intersection(con_b)))
					# cc_v = float(len(con_b.intersection(con_a)))

					fcc, fcc_v = cc * 1.0/len(con_a), cc * 1.0/len(con_b)
					self.fcc_matrix.append(f'{i + 1} {j + 1} {fcc:.3f} {fcc_v:.3f}')

			with open(self.fcc_matrix_f,'w') as out:
				out.write('\n'.join(self.fcc_matrix))
			out.close()
		elif not self.fcc_matrix:

			# print('+ Loading FCC matrix')

			with open(self.fcc_matrix_f) as f:
				for l in f.readlines():
					self.fcc_matrix.append(l)
			f.close()
		else:

			# print('+ Using previously loaded FCC matrix')

			pass

	# def calculate_rmsd(self):
	# 	""" Calculate RMSD using lowest HADDOCK-score structure as reference """
	#
	# 	hs_list = [(k, self.structure_dic[k]['haddock-score']) for k in self.structure_dic]
	# 	sorted_hs_list = sorted(hs_list, key=lambda x: (-x[1], x[0]))
	# 	sorted_hs_list.reverse()
	# 	sorted_pdb_list = [e[0] for e in sorted_hs_list]
	#
	# 	rmsd_dic = profit_rmsd(sorted_pdb_list)

	def run_dockq(self):

		print('+ Running DockQ')

		if self.param_dic['input']['reference'] == 'lowest':
			score_list = []
			for pdb in self.structure_dic:
				score_list.append((pdb, self.structure_dic[pdb]['haddock-score']))
			sorted_score_list = sorted(score_list, key=lambda x: x[1])
			reference_pdb = sorted_score_list[0][0]
			reference_score = sorted_score_list[0][1]

			print(f'++ Using {reference_pdb} as reference structure, lowest haddock score: {reference_score:.3f}')

		else:
			reference_pdb = self.param_dic['input']['reference']

			print(f'++ Using {reference_pdb} as reference structure, user input')

		for pdb in self.structure_dic:
			result_dic = dockq(reference_pdb, pdb)
			self.structure_dic[pdb]['dockq'] = result_dic


class Cluster(object):
	"""Defines a Cluster. A Cluster is created with a name and a center (Element class)"""

	__slots__ = ['name', 'center', 'members']

	def __init__(self, name, center):
		self.name = name
		self.center = center

		self.members = []

		self.populate()

	def __len__(self):
		return len(self.members) + 1  # +1 Center

	def populate(self):
		"""
		Populates the Cluster member list through the
		neighbor list of its center.
		"""

		name = self.name
		# Assign center
		ctr = self.center
		ctr.assign_cluster(name)

		mlist = self.members
		# Assign members
		ctr_nlist = (n for n in ctr.neighbors if not n.cluster)
		for e in ctr_nlist:
			mlist.append(e)
			e.assign_cluster(name)

	def add_member(self, element):
		"""
		Adds one single element to the cluster.
		"""
		l = self.members
		l.append(element)
		element.assign_cluster(self.name)


class Element(object):
	"""Defines a 'clusterable' Element"""

	__slots__ = ['name', 'cluster', 'neighbors']

	def __init__(self, name):
		self.name = name
		self.cluster = 0
		self.neighbors = set()

	def add_neighbor(self, neighbor):
		"""Adds another element to the neighbor list"""
		self.neighbors.add(neighbor)

	def assign_cluster(self, clust_id):
		"""Assigns the Element to Cluster. 0 if unclustered"""
		self.cluster = clust_id
