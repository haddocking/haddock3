import glob
# import operator
import re
import os

from haddock.workflows.scoring.analysis.contact import run_contacts
from haddock.workflows.scoring.analysis.dfire import dfire
from haddock.workflows.scoring.analysis.fastcontact import fastcontact
# from haddock.workflows.scoring.src import calc_fcc_matrix, cluster_fcc
from haddock.workflows.scoring.analysis.profit import profit_rmsd
from haddock.workflows.scoring.src import calc_fcc_matrix


class Ana:

	def __init__(self):
		self.structure_dic = {}
		self.contact_filelist = []
		if not os.path.isdir('contacts'):
			os.system('mkdir contacts')

		self.fcc_matrix_f = 'fcc.matrix'

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
		for pdb in self.structure_dic:
			fast_elec, fast_desol = fastcontact(pdb)
			self.structure_dic[pdb]['fastelec'] = fast_elec
			self.structure_dic[pdb]['fastdesol'] = fast_desol

	def run_dfire(self):
		""" Run dfire on all scored PDBs """
		for pdb in self.structure_dic:
			d_binding, d_score = dfire(pdb)
			self.structure_dic[pdb]['dfire-ebinding'] = d_binding
			self.structure_dic[pdb]['dfire-score'] = d_score

	def calculate_contacts(self):
		""" Calculate atomic contacts """

		for pdb in self.structure_dic:
			contact_out = run_contacts(pdb)
			self.contact_filelist.append(contact_out)

	# TODO: Implement FCC clustering
	def calc_fcc_matrix(self):
		pass

	def calculate_rmsd(self):
		""" Calculate RMSD using lowest HADDOCK-score structure as reference """

		hs_list = [(k, self.structure_dic[k]['haddock-score']) for k in self.structure_dic]
		sorted_hs_list = sorted(hs_list, key=lambda x: (-x[1], x[0]))
		sorted_hs_list.reverse()
		sorted_pdb_list = [e[0] for e in sorted_hs_list]

		rmsd_dic = profit_rmsd(sorted_pdb_list)


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
