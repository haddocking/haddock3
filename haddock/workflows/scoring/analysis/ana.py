# import glob
import re
import os
import haddock.workflows.scoring.config as config
from haddock.workflows.scoring.analysis.contact import Contacts
# from haddock.workflows.scoring.analysis.dfire import dfire
# from haddock.workflows.scoring.analysis.dockq import dockq
# from haddock.workflows.scoring.analysis.fastcontact import fastcontact


class Ana:

	def __init__(self, pdb_l):
		self.structure_dic = {}
		self.contact_filelist = []
		self.structure_haddockscore_list = []
		self.fcc_matrix_f = 'fcc.matrix'
		self.output_f = 'scoring.stat'
		self.fcc_matrix = []
		self.con_list = None

		for pdb in pdb_l:
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
					# dirty account for 8.754077E-02 notation and the eventual $DANI or $NOE
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
					# WARNING: This is highly dependent on the values printed by the recipe, refer to print_coorheader.cns
					total, bonds, angles, improper, dihe, vdw, elec, air, cdih, coup, rdcs, vean, dani, xpcs, rg = energy_v
					vdw = float(vdw)
					elec = float(elec)
					air = float(air)

				if 'REMARK Desolvation' in line:
					desolv = float(re.findall(desolv_regex, line)[0])

				if 'REMARK buried surface area' in line:
					bsa = float(re.findall(bsa_regex, line)[0])
					break

			f.close()

			self.structure_dic[pdb]['haddock-score'] = None
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

		try:
			for pdb in self.structure_dic:
				_ = self.structure_dic[pdb]['vdw']
		except KeyError:
			self.extract_energies()

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

		wd = os.getcwd()
		output_f = f'{wd}/cluster_{cutoff}_{threshold}.out'

		self.calculate_contacts()
		self.calc_fcc_matrix()

		cutoff = float(cutoff)
		partner_cutoff = float(cutoff) * float(strictness)

		# populate internal data structure with cluster info
		for pdb in self.structure_dic:
			self.structure_dic[pdb][f"cluster-{cutoff}-{threshold}_name"] = float('nan')
			self.structure_dic[pdb][f"cluster-{cutoff}-{threshold}_internal_ranking"] = float('nan')
			self.structure_dic[pdb][f"cluster-{cutoff}-{threshold}_overall_ranking"] = float('nan')

		elements = {}

		for e in self.fcc_matrix:
			ref, mobi, d_rm, d_mr = e.split()
			ref = int(ref)
			mobi = int(mobi)
			d_rm = float(d_rm)
			d_mr = float(d_mr)

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
			if d_rm >= cutoff and d_mr >= partner_cutoff:
				r.add_neighbor(m)
			if d_mr >= cutoff and d_rm >= partner_cutoff:
				m.add_neighbor(r)

		clusters = []
		# threshold -= 1  # Account for center
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

		path = '/'.join(list(self.structure_dic.items())[0][0].split('/')[:-1])
		cluster_dic = {}
		for c in clusters:
			clustered_models_list = []

			# add cluster center to internal structure
			# cluster_center_str = '0' * (6 - len(str(c.center.name - 1))) + str(c.center.name - 1)
			cluster_center_str = '0' * (6 - len(str(c.center.name))) + str(c.center.name)
			cluster_center_name = f'{path}/{cluster_center_str}_conv.pdb'
			self.structure_dic[cluster_center_name][f"cluster-{cutoff}-{threshold}_name"] = c.name
			self.structure_dic[cluster_center_name][f"cluster-{cutoff}-{threshold}_internal_ranking"] = 0

			sort_tag = True
			for model in [m.name for m in c.members]:
				# model_str = '0' * (6 - len(str(model - 1))) + str(model - 1)
				model_str = '0' * (6 - len(str(model))) + str(model)
				name = f'{path}/{model_str}_conv.pdb'
				try:
					haddock_score = self.structure_dic[name]['haddock-score']
				except KeyError:
					# No haddock score...
					haddock_score = float('nan')
					sort_tag = False

				clustered_models_list.append((model, haddock_score))

				# add cluster info to internal structure
				self.structure_dic[name][f"cluster-{cutoff}-{threshold}_name"] = c.name
			if sort_tag:
				# sort cluster elements by haddock score
				model_list = sorted(clustered_models_list, key=lambda x: x[1])
			else:
				# sort by model number
				model_list = sorted(clustered_models_list, key=lambda x: x[0])
			score_list = [e[1] for e in clustered_models_list]
			mean_score = sum(score_list) / len(score_list)
			# top4_mean_score = sum(score_list[:4]) / len(score_list[:4])
			tbw = f"Cluster {c.name} -> ({c.center.name}) "
			for i, m in enumerate(model_list):
				tbw += f'{m[0]} '
				# model_str = '0' * (6 - len(str(m[0] - 1))) + str(m[0] - 1)
				model_str = '0' * (6 - len(str(m[0]))) + str(m[0])
				name = f'{path}/{model_str}_conv.pdb'
				self.structure_dic[name][f"cluster-{cutoff}-{threshold}_internal_ranking"] = i

			tbw += '\n'
			cluster_dic[c.name] = (tbw, mean_score)

		# sort clusters by mean haddock score
		sorted_cluster_list = sorted([(x, cluster_dic[x][1]) for x in list(cluster_dic.keys())], key=lambda x: x[1])

		cluster_ranking = dict([(e[0], i + 1) for i, e in enumerate(sorted_cluster_list)])
		# add overall cluster ranking to data structure
		for pdb in self.structure_dic:
			cluster_name = self.structure_dic[pdb][f"cluster-{cutoff}-{threshold}_name"]
			if cluster_name == cluster_name:  # means it is not NaN
				overall_ranking = cluster_ranking[cluster_name]
				self.structure_dic[pdb][f"cluster-{cutoff}-{threshold}_overall_ranking"] = overall_ranking

		# output
		with open(output_f, 'w') as out:
			for c in sorted_cluster_list:
				cluster_name, cluster_mean = c
				tbw_l = cluster_dic[cluster_name][0].split('->')
				tbw = f'{tbw_l[0]}[{cluster_mean:.3f}] ->{tbw_l[1]}'
				out.write(tbw)
		out.close()

		return output_f

	def calculate_contacts(self):
		""" Calculate atomic contacts """

		pdb_list = list(self.structure_dic.keys())

		con = Contacts()
		con.calculate_contacts(pdb_list, cutoff=5.0)

		self.con_list = con.contact_file_list

	def calc_fcc_matrix(self):
		""" Calculate the FCC matrix (extracted and adapted from calc_fcc_matrix.py """

		# if not os.path.isfile(self.fcc_matrix_f) and not self.fcc_matrix:

		# print('+ Creating FCC matrix')
		if not self.con_list:
			self.calculate_contacts()

		# very important!
		self.con_list.sort()

		contacts = []
		for con_f in self.con_list:
			t = []
			with open(con_f) as f:
				for l in f.readlines():
					t.append(int(l))
			f.close()
			contacts.append(set(t))

		# get the pairwise matrix
		for i in range(len(contacts)):
			for j in range(i + 1, len(contacts)):
				con_a = contacts[i]
				con_b = contacts[j]

				cc = float(len(con_a.intersection(con_b)))
				# cc_v = float(len(con_b.intersection(con_a)))

				fcc, fcc_v = cc * 1.0 / len(con_a), cc * 1.0 / len(con_b)
				self.fcc_matrix.append(f'{i + 1} {j + 1} {fcc:.3f} {fcc_v:.3f}')
		#
		# 	with open(self.fcc_matrix_f, 'w') as out:
		# 		out.write('\n'.join(self.fcc_matrix))
		# 	out.close()
		# elif not self.fcc_matrix:
		#
		# 	# print('+ Loading FCC matrix')
		#
		# 	with open(self.fcc_matrix_f) as f:
		# 		for l in f.readlines():
		# 			self.fcc_matrix.append(l)
		# 	f.close()
		# else:
		#
		# 	# print('+ Using previously loaded FCC matrix')
		# 	pass

	def run_dockq(self):

		print('+ Running DockQ')

		if config.param_dic['input']['reference'] == 'lowest':
			score_list = []
			for pdb in self.structure_dic:
				score_list.append((pdb, self.structure_dic[pdb]['haddock-score']))
			sorted_score_list = sorted(score_list, key=lambda x: x[1])
			reference_pdb = sorted_score_list[0][0]
			reference_score = sorted_score_list[0][1]

			print(f'++ Using {reference_pdb} as reference structure, lowest haddock score: {reference_score:.3f}')

		else:
			reference_pdb = config.param_dic['input']['reference']

			print(f'++ Using {reference_pdb} as reference structure, user input')

		# pdb_list = list(self.structure_dic.keys())
		#
		# dockq = DockQ()
		# d = dockq.run(pdb_list, reference_pdb)

		#
		for pdb in self.structure_dic:
			result_dic = dockq(reference_pdb, pdb)
			for k in result_dic:
				self.structure_dic[pdb][k] = result_dic[k]

	def output(self):

		print(f'+ Saving results to {self.output_f}')

		# sort by haddock score!
		score_list = [(pdb, self.structure_dic[pdb]['haddock-score']) for pdb in self.structure_dic]
		sorted_score_list = sorted(score_list, key=lambda x: x[1])
		k = sorted_score_list[0][0]
		header = 'model ranking ' + ' '.join(list(self.structure_dic[k])) + '\n'

		with open(self.output_f, 'w') as out:
			out.write(header)

			for i, e in enumerate(sorted_score_list):
				pdb = e[0]
				tbw = f'{pdb} {i+1} '
				for v in self.structure_dic[pdb].values():
					if isinstance(v, float):
						tbw += f'{v:.3f} '
					else:
						tbw += f'{v} '
				tbw += '\n'
				out.write(tbw)
		out.close()


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
