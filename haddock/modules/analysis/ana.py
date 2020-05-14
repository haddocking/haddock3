# import glob
import re
import os
import shutil
import statistics
import math
import subprocess
import itertools
import configparser
import operator
from haddock.modules.structure.utils import PDB
from haddock.utils.files import get_full_path
from haddock.modules.analysis.fcc import read_matrix, cluster_elements, output_clusters
from haddock.modules.analysis.contact import Contacts
from haddock.modules.analysis.dfire import dfire
from haddock.modules.analysis.dockq import dockq
from haddock.modules.analysis.fastcontact import fastcontact

etc_folder = get_full_path('haddock', 'etc')
config_file = os.path.join(etc_folder, 'haddock3.ini')

ini = configparser.ConfigParser(os.environ)
ini.read(config_file, encoding='utf-8')


class Ana:

	def __init__(self, pdb_l):
		self.structure_dic = {}
		self.contact_filelist = []
		self.structure_haddockscore_list = []
		self.fcc_matrix_f = 'fcc.matrix'
		self.ss_output_f = 'ss'
		self.cluster_output_f = 'cluster'
		self.fcc_matrix = []
		self.con_list = None

		for pdb in pdb_l:
			self.structure_dic[pdb] = {}

		if not os.path.isdir('analysis'):
			os.mkdir('analysis')

	def extract_energies(self):
		""" Extract energies from the header of the PDB file, according to HADDOCK formatting """

		total = .0
		bonds = .0
		angles = .0
		improper = .0
		dihe = .0
		vdw = .0
		elec = .0
		air = .0
		cdih = .0
		coup = .0
		sani = .0
		vean = .0
		dani = .0
		desolv = .0
		bsa = .0

		for pdb in self.structure_dic:

			vdw_elec_air_regex = r"\s(\-?\d*\.?\d{1,}(E-\d{1,}|)|0\b|\$\w{1,})"
			desolv_regex = r"(\-?\d*\.?\d*)$"
			bsa_regex = r"(\-?\d*\.?\d*)$"
			inter_regex = r"(\-?\d*\.?\d*)$"

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

				# if 'REMARK Internal' in line:
				# 	internal = float(re.findall(inter_regex, line)[0])

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
			# self.structure_dic[pdb]['internal']

	def calculate_haddock_score(self):
		""" Calculate the HADDOCK Score of the PDB file using its appropriate weight """

		print(f'\n+ Calculating HADDOCK score for {len(self.structure_dic)} structures')

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

	def cluster(self, cutoff, size=4, strictness=0.75, detailed=False):
		""" Cluster scored models using FCC, output a sorted text file containing clusters and mean scores """

		print(f'+ FCC clustering with cutoff: {cutoff:.2f} and max size: {size}')

		self.calculate_contacts()  # this will be sorted by HADDOCK Score
		fcc_matrix = self.calc_fcc_matrix()
		pool = read_matrix(fcc_matrix, cutoff, strictness)
		element_pool, clusters = cluster_elements(pool, size)
		# output_clusters(f'cluster_debug_{cutoff:.2f}-{int(size)}.out', clusters)

		# populate internal data structure with cluster info
		for pdb in self.structure_dic:
			self.structure_dic[pdb][f"cluster-{cutoff}-{size}_name"] = .0
			self.structure_dic[pdb][f"cluster-{cutoff}-{size}_internal_ranking"] = .0
			self.structure_dic[pdb][f"cluster-{cutoff}-{size}_overall_ranking"] = .0

		# Important: guess the root, this might change according to the worflow
		first_model_name = list(self.structure_dic.items())[0][0]
		path_regex = r"(.*)/"
		root_regex = r"/(.*)_\d*\.pdb"

		path = re.findall(path_regex, first_model_name)[0]
		root = re.findall(root_regex, first_model_name)[0]
		wd = os.getcwd()
		output_f = f'{wd}/analysis/fcc_{cutoff}-{size}.out'
		# output_f = f'{wd}/{path}/analysis/fcc_{cutoff}-{threshold}.out'

		# get the ORDER, +1 since line counting starts at 1
		ref_dic = dict([(i+1, int(e.split('/')[-1].split('.')[0].split('_')[-1])) for i, e in enumerate(self.con_list)])

		cluster_dic = {}
		for c in clusters:
			# The numbers are the INDEX of the INPUT FILE, not the MODEL numbers
			c.center.name -= 1
			clustered_models_list = []

			center_model_number = ref_dic[c.center.name]

			cluster_center_name = f'{path}/{root}_{center_model_number}.pdb'
			# cluster_center_name = f'{path}/{root}_{c.center.name}.pdb'
			if not os.path.isfile(cluster_center_name):
				cluster_center_name_wzeros = f'{path}/{root}_{str(center_model_number).zfill(6)}.pdb'
				if os.path.isfile(cluster_center_name_wzeros):
					cluster_center_name = cluster_center_name_wzeros

			self.structure_dic[cluster_center_name][f"cluster-{cutoff}-{size}_name"] = c.name
			self.structure_dic[cluster_center_name][f"cluster-{cutoff}-{size}_internal_ranking"] = 0

			sort_tag = True
			for model_index in [m.name for m in c.members]:
				model_number = ref_dic[model_index]
				name = f'{path}/{root}_{model_number}.pdb'
				if not os.path.isfile(name):
					name_wzeros = f'{path}/{root}_{str(model_number).zfill(6)}.pdb'
					if os.path.isfile(name_wzeros):
						name = name_wzeros

				haddock_score = .0
				try:
					haddock_score = self.structure_dic[name]['haddock-score']
				except KeyError:
					# No haddock score...?
					print('+ ERROR: Model has no Haddock Score')
					exit()
				# haddock_score = .0
				# sort_tag = False

				clustered_models_list.append((model_number, haddock_score))

				# add cluster info to internal structure
				self.structure_dic[name][f"cluster-{cutoff}-{size}_name"] = c.name

			if sort_tag:
				# sort cluster elements by haddock score
				model_list = sorted(clustered_models_list, key=lambda x: x[1])
			else:
				# sort by model number
				model_list = sorted(clustered_models_list, key=lambda x: x[0])

			score_list = [e[1] for e in clustered_models_list]
			mean_score = sum(score_list) / len(score_list)
			top4_mean_score = sum(score_list[:4]) / len(score_list[:4])
			cluster_dic[c.name] = (mean_score, top4_mean_score, c.center.name, [])

			for i, m in enumerate(model_list):
				model_id = m[0]

				cluster_dic[c.name][3].append(model_id)
				name = f'{path}/{root}_{model_id}.pdb'

				if not os.path.isfile(name):
					name_wzeros = f'{path}/{root}_{str(model_id).zfill(6)}.pdb'
					if os.path.isfile(name_wzeros):
						name = name_wzeros

				self.structure_dic[name][f"cluster-{cutoff}-{size}_internal_ranking"] = i+1

		# sort clusters by mean haddock score
		sorted_cluster_list = sorted(
				[(x, cluster_dic[x][0], cluster_dic[x][1], cluster_dic[x][2], cluster_dic[x][3]) for x in list(cluster_dic.keys())],
				key=lambda x: x[1])

		# sort cluster by top4 score
		# sorted_cluster_list = sorted([(x, cluster_dic[x][0], cluster_dic[x][1], cluster_dic[x][2], cluster_dic[x][3) for x in list(cluster_dic.keys())], key=lambda x: x[2])

		cluster_ranking = dict([(e[0], i + 1) for i, e in enumerate(sorted_cluster_list)])
		# add overall cluster ranking to data structure
		for pdb in self.structure_dic:
			cluster_name = self.structure_dic[pdb][f"cluster-{cutoff}-{size}_name"]
			if cluster_name != 0:
				overall_ranking = cluster_ranking[cluster_name]
				self.structure_dic[pdb][f"cluster-{cutoff}-{size}_overall_ranking"] = overall_ranking

		# output
		# header = '# Cluster ID\t[mean_score, top4_mean_score]\t->\t(center) top1 top2 top3 top4 ...\n'
		with open(output_f, 'w') as out:
			# out.write(header)
			for c in sorted_cluster_list:
				cluster_name, cluster_mean, cluster_top4_mean, cluster_center, sorted_model_list = c
				sorted_model_string = ' '.join(map(str, sorted_model_list))
				tbw = f'Cluster {cluster_name} [{cluster_mean:.2f}, {cluster_top4_mean:.2f}] -> ({cluster_center}) ' \
					f'{sorted_model_string}\n'
				# tbw = f'Cluster {cluster_name} [{cluster_mean:.2f}, {cluster_top4_mean:.2f}] -> {sorted_model_string}\n'
				out.write(tbw)
		out.close()

		if detailed:
			# write detailed output files
			total_cluster = len(sorted_cluster_list)

			dh_dic = {j + 1: [] for j in range(total_cluster)}
			edesol_dic = {j + 1: [] for j in range(total_cluster)}

			component_dic = {'size': None, 'inter': [], 'nb': [], 'vdw': [], 'elec': [], 'air': [], 'cdih': [], 'coup': [], 'sani': [], 'vean': [], 'dani': []}
			ener_dic = {j + 1: component_dic for j in range(total_cluster)}

			haddock_score_dic = {j + 1: [] for j in range(total_cluster)}
			bsa_dic = {j + 1: [] for j in range(total_cluster)}

			rmsd_dic = {j + 1: [] for j in range(total_cluster)}
			rmsdmin_dic = {j + 1: [] for j in range(total_cluster)}
			viol_dic = {j + 1: [] for j in range(total_cluster)}

			stat_dic = {'haddock-score': [], 'rmsd': [], 'rmsd-Emin': [], 'size': None, 'inter': [], 'nb': [], 'vdw': [], 'elec': [], 'air': [], 'cdih': [], 'coup': [], 'sani': [], 'vean': [], 'dani': []}

			cluster_stat_dic = {j + 1: stat_dic for j in range(total_cluster)}
			cluster_stat_best2_dic = {j + 1: stat_dic for j in range(total_cluster)}
			cluster_stat_best4_dic = {j + 1: stat_dic for j in range(total_cluster)}

			for c in sorted_cluster_list:
				cluster_name, cluster_mean, cluster_top4_mean, cluster_center, sorted_model_list = c

				file_cns_tbw = f''
				file_cns_tbw_best4 = f''
				file_cns_tbw_best2 = f''
				file_nam_tbw = f''
				file_nam_tbw_best2 = f''
				file_nam_tbw_best4 = f''
				file_list_tbw = f''
				file_list_tbw_best2 = f''
				file_list_tbw_best4 = f''

				file_nam_bsa = f'#struc BSA\n'
				file_nam_dh = f'#struc dH\n'
				file_nam_desol = f'#struc Edesolv\n'
				# file_nam_ener = f'#struc Einter Enb Evdw+0.1Eelec Evdw Eelec Eair Ecdih Ecoup Esani Evean Edani\n' # old
				file_nam_ener = f'#struc Einter Evdw Eelec Eair Ecdih Ecoup Esani Evean Edani\n'  # new
				file_nam_haddockscore = f'#struc haddock-score\n'
				file_nam_rmsd = f'#struc rmsd_all\n'
				file_nam_rmsdmin = f'#struc rmsd_Emin\n'
				file_nam_viol = f'#struc #NOEviol #Dihedviol #Coupviol #Saniviol #Veanviol #Daniviol\n'

				ener_dic[cluster_name]['size'] = len(sorted_model_list)

				for i, model_id in enumerate(sorted_model_list):
					name = f'{path}/{root}_{model_id}.pdb'

					haddock_score_v = self.structure_dic[name]['haddock-score']
					bsa_v = self.structure_dic[name]['bsa']
					desolv_v = self.structure_dic[name]['desolv']
					inter_v = self.structure_dic[name]['total']
					vdw_v = self.structure_dic[name]['vdw']
					elec_v = self.structure_dic[name]['elec']
					air_v = self.structure_dic[name]['air']
					cdih_v = self.structure_dic[name]['cdih']
					coup_v = self.structure_dic[name]['coup']
					sani_v = self.structure_dic[name]['sani']
					vean_v = self.structure_dic[name]['vean']
					dani_v = self.structure_dic[name]['dani']

					# cluster_bsa.txt
					bsa_dic[cluster_name].append(bsa_v)

					# PENDING #============================================================================================#
					# cluster_dH.txt
					# dh_dic[cluster_name].append(dh)
					# =====================================================================================================#

					# cluster_Edesolv.txt
					edesol_dic[cluster_name].append(desolv_v)

					# cluster_ener.txt - contains the average energy terms of each cluster and the standard deviations
					ener_dic[cluster_name]['inter'].append(inter_v)
					ener_dic[cluster_name]['vdw'].append(vdw_v)
					ener_dic[cluster_name]['elec'].append(elec_v)
					ener_dic[cluster_name]['air'].append(air_v)
					ener_dic[cluster_name]['cdih'].append(cdih_v)
					ener_dic[cluster_name]['coup'].append(coup_v)
					ener_dic[cluster_name]['sani'].append(sani_v)
					ener_dic[cluster_name]['vean'].append(vean_v)
					ener_dic[cluster_name]['dani'].append(dani_v)

					# cluster_haddock.txt - contains the average combined haddock score
					haddock_score_dic[cluster_name].append(haddock_score_v)

					# PENDING #============================================================================================#
					# cluster_rmsd.txt - contains the average RMSD and standard deviation from the best (lowest)
					#   HADDOCK score structure of cluster of the structures belonging to that cluster
					# rmsd_dic[cluster_name].append(rmsd)
					#
					# cluster_rmsd-Emin.txt - contains the average RMSD and standard deviation of the clusters from the best
					#   (lowest) HADDOCK score structure of all calculated structures
					# rmsdmin_dic[cluster_name].append(rmsdmin)
					#
					# cluster_viol.txt -  contains the average AIR and dihedral violations for each cluster and the
					#   standard deviations
					# viol_dic[cluster_name].append(viol)
					#
					# file.nam_clustX_dH - contains the total energy difference calculated as total energy of the complex -
					#   Sum of total energies of the individual components
					# file_nam_dh += ''
					#
					# file.nam_clustX_rmsd - contains the RMSD of each structure of cluster X from the best (lowest) HADDOCK
					#   score structure of cluster X.
					# file_nam_rmsd += f''
					#
					# file.nam_clustX_rmsd-Emin - contains the RMSD of each structure of cluster X from the best (lowest)
					#   HADDOCK score structure of all calculated structures
					# file_nam_rmsdmin += f''
					#
					# file.nam_clustX_viol - contains the number of AIR and dihedral violations per structure
					# Note: The ordering of the structures in those files follows the HADDOCK score ranking.
					# file_nam_viol += f''
					# =====================================================================================================#

					# DEPRECATED #=========================================================================================#
					# # file.cns_clustX - contains the name of all the pdb files that belong to the cluster X (CNS format)
					# file_cns_tbw += f'evaluate (&filenames.bestfile_{i}="{name}")\n'
					# # file.cns_clustX_bestY - contains the name of the best Y pdb files that belong to the cluster X
					# #   (CNS format)
					# if i < 2:
					# 	file_cns_tbw_best2 += f'evaluate (&filenames.bestfile_{i}="{name}")\n'
					# if i < 4:
					# 	file_cns_tbw_best4 += f'evaluate (&filenames.bestfile_{i}="{name}")\n'
					#
					# # file.nam_clustX - contains the name of all the pdb files that belong to the cluster X
					# file_nam_tbw += f'{name}\n'
					# # file.nam_clustX_bestY - contains the name of the best Y pdb files that belong to the cluster X
					# if i < 2:
					# 	file_nam_tbw_best2 += f'{name}\n'
					# if i < 4:
					# 	file_nam_tbw_best4 += f'{name}\n'
					# =====================================================================================================#

					# file.list_clustX
					file_list_tbw += f'{name} {{ {haddock_score_v:.2f} }}\n'
					# file.list_clustX_bestY
					if i < 2:
						file_list_tbw_best2 += f'{name} {{ {haddock_score_v:.2f} }}\n'
					if i < 4:
						file_list_tbw_best4 += f'{name} {{ {haddock_score_v:.2f} }}\n'

					# file.nam_clustX_bsa
					file_nam_bsa += f'{name} {bsa_v:.2f}\n'

					# file.nam_clustX_Edesol
					file_nam_desol += f'{name} {desolv_v:.2f}\n'

					# file.nam_clustX_ener
					file_nam_ener += f'{name} {inter_v:.2f} {vdw_v:.2f} {elec_v:.2f} {air_v:.2f} {cdih_v:.2f} {sani_v:.2f} {vean_v:.2f} {dani_v:.2f}\n'

					# file.nam_clustX_haddock-score
					file_nam_haddockscore += f'{name} {haddock_score_v:.2f}\n'

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}_bsa', 'w') as f:
					f.write(file_nam_bsa)
				f.close()

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}_dH', 'w') as f:
					f.write(file_nam_dh)
				f.close()

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}_Edesolv', 'w') as f:
					f.write(file_nam_desol)
				f.close()

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}_ener', 'w') as f:
					f.write(file_cns_tbw)
				f.close()

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}_haddock-score', 'w') as f:
					f.write(file_nam_haddockscore)
				f.close()

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}', 'w') as f:
					f.write(file_list_tbw)
				f.close()

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}_best2', 'w') as f:
					f.write(file_list_tbw_best2)
				f.close()

				with open(f'{wd}/analysis/fcc_{cutoff}-{size}_clust{cluster_name}_best4', 'w') as f:
					f.write(file_list_tbw_best4)
				f.close()

				# PENDING #================================================================================================#
				# with open(f'{wd}/analysis/file.nam_clust{cluster_name}_rmsd', 'w') as f:
				# 	f.write(file_nam_rmsd)
				# f.close()
				#
				# with open(f'{wd}/analysis/file.nam_clust{cluster_name}_rmsd-Emin', 'w') as f:
				# 	f.write(file_nam_rmsdmin)
				# f.close()
				#
				# with open(f'{wd}/analysis/file.cns_clust{cluster_name}_viol', 'w') as f:
				# 	f.write(file_nam_viol)
				# f.close()
				# =========================================================================================================#

				# DEPRECATED #=============================================================================================#
				# with open(f'{wd}/analysis/file.cns_clust{cluster_name}', 'w') as f:
				# 	f.write(file_cns_tbw)
				# f.close()
				#
				# with open(f'{wd}/analysis/file.cns_clust{cluster_name}_best2', 'w') as f:
				# 	f.write(file_cns_tbw_best2)
				# f.close()
				#
				# with open(f'{wd}/analysis/file.cns_clust{cluster_name}_best4', 'w') as f:
				# 	f.write(file_cns_tbw_best4)
				# f.close()
				#
				# with open(f'{wd}/analysis/file.nam_clust{cluster_name}', 'w') as f:
				# 	f.write(file_nam_tbw)
				# f.close()
				#
				# with open(f'{wd}/analysis/file.nam_clust{cluster_name}_best2', 'w') as f:
				# 	f.write(file_nam_tbw_best2)
				# f.close()
				#
				# with open(f'{wd}/analysis/file.nam_clust{cluster_name}_best4', 'w') as f:
				# 	f.write(file_nam_tbw_best4)
				# f.close()
				# =========================================================================================================#

			with open(f'{wd}/analysis/fcc_{cutoff}-{size}_bsa.txt', 'w') as f:
				f.write('#Cluster BSA sd\n')
				for cluster_id in bsa_dic:
					f.write(f'file.nam_clust{cluster_id} {statistics.mean(bsa_dic[cluster_id]):.2f} {statistics.stdev(bsa_dic[cluster_id]):.2f}\n')
			f.close()

			with open(f'{wd}/analysis/fcc_{cutoff}-{size}_haddock.txt', 'w') as f:
				f.write('#Cluster haddock-score sd\n')
				for cluster_id in haddock_score_dic:
					f.write(f'file.nam_clust{cluster_id} {statistics.mean(haddock_score_dic[cluster_id]):.2f} {statistics.stdev(haddock_score_dic[cluster_id]):.2f}\n')
			f.close()

			with open(f'{wd}/analysis/fcc_{cutoff}-{size}.stats', 'w') as f:
				f.write('#Cluster haddock-score sd rmsd sd rmsd-Emin sd Nstruc Etot sd Evdw sd Eelec sd BSA sd\n')
				for cluster_id in range(1, total_cluster+1):

					hs_mean = statistics.mean(haddock_score_dic[cluster_id])
					hs_sd = statistics.stdev(haddock_score_dic[cluster_id])

					rmsd_mean = .0
					rmsd_sd = .0

					rmsdemin_mean = .0
					rmsdemin_sd = .0

					size = len(haddock_score_dic[cluster_id])

					total_mean = statistics.mean(ener_dic[cluster_id]['inter'])
					total_sd = statistics.stdev(ener_dic[cluster_id]['inter'])

					vdw_mean = statistics.mean(ener_dic[cluster_id]['vdw'])
					vdw_sd = statistics.stdev(ener_dic[cluster_id]['vdw'])

					elec_mean = statistics.mean(ener_dic[cluster_id]['elec'])
					elec_sd = statistics.stdev(ener_dic[cluster_id]['elec'])

					bsa_mean = statistics.mean(bsa_dic[cluster_id])
					bsa_sd = statistics.stdev(bsa_dic[cluster_id])

					f.write(f'file.nam_clust{cluster_id} {hs_mean:.2f} {hs_sd:.2f} {rmsd_mean:.2f} {rmsd_sd:.2f} {rmsdemin_mean:.2f} {rmsdemin_sd:.2f} {size} {total_mean:.2f} {total_sd:.2f} {vdw_mean:.2f} {vdw_sd:.2f} {elec_mean:.2f} {elec_sd:.2f} {bsa_mean:.2f} {bsa_sd:.2f}\n')
			f.close()

			with open(f'{wd}/analysis/fcc_{cutoff}-{size}.stats_best2', 'w') as f:
				f.write('#Cluster haddock-score sd rmsd sd rmsd-Emin sd Nstruc Etot sd Evdw sd Eelec sd BSA sd\n')
				for cluster_id in range(1, total_cluster + 1):

					hs_mean = statistics.mean(haddock_score_dic[cluster_id][:2])
					hs_sd = statistics.stdev(haddock_score_dic[cluster_id][:2])

					rmsd_mean = .0
					rmsd_sd = .0

					rmsdemin_mean = .0
					rmsdemin_sd = .0

					size = len(haddock_score_dic[cluster_id])

					total_mean = statistics.mean(ener_dic[cluster_id]['inter'][:2])
					total_sd = statistics.stdev(ener_dic[cluster_id]['inter'][:2])

					vdw_mean = statistics.mean(ener_dic[cluster_id]['vdw'][:2])
					vdw_sd = statistics.stdev(ener_dic[cluster_id]['vdw'][:2])

					elec_mean = statistics.mean(ener_dic[cluster_id]['elec'][:2])
					elec_sd = statistics.stdev(ener_dic[cluster_id]['elec'][:2])

					bsa_mean = statistics.mean(bsa_dic[cluster_id][:2])
					bsa_sd = statistics.stdev(bsa_dic[cluster_id][:2])

					f.write(f'file.nam_clust{cluster_id} {hs_mean:.2f} {hs_sd:.2f} {rmsd_mean:.2f} {rmsd_sd:.2f} {rmsdemin_mean:.2f} {rmsdemin_sd:.2f} {size} {total_mean:.2f} {total_sd:.2f} {vdw_mean:.2f} {vdw_sd:.2f} {elec_mean:.2f} {elec_sd:.2f} {bsa_mean:.2f} {bsa_sd:.2f}\n')
			f.close()

			with open(f'{wd}/analysis/fcc_{cutoff}-{size}.stats_best4', 'w') as f:
				f.write('#Cluster haddock-score sd rmsd sd rmsd-Emin sd Nstruc Etot sd Evdw sd Eelec sd BSA sd\n')
				for cluster_id in range(1, total_cluster + 1):

					hs_mean = statistics.mean(haddock_score_dic[cluster_id][:4])
					hs_sd = statistics.stdev(haddock_score_dic[cluster_id][:4])

					rmsd_mean = .0
					rmsd_sd = .0

					rmsdemin_mean = .0
					rmsdemin_sd = .0

					size = len(haddock_score_dic[cluster_id])

					total_mean = statistics.mean(ener_dic[cluster_id]['inter'][:4])
					total_sd = statistics.stdev(ener_dic[cluster_id]['inter'][:4])

					vdw_mean = statistics.mean(ener_dic[cluster_id]['vdw'][:4])
					vdw_sd = statistics.stdev(ener_dic[cluster_id]['vdw'][:4])

					elec_mean = statistics.mean(ener_dic[cluster_id]['elec'][:4])
					elec_sd = statistics.stdev(ener_dic[cluster_id]['elec'][:4])

					bsa_mean = statistics.mean(bsa_dic[cluster_id][:4])
					bsa_sd = statistics.stdev(bsa_dic[cluster_id][:4])

					f.write(f'file.nam_clust{cluster_id} {hs_mean:.2f} {hs_sd:.2f} {rmsd_mean:.2f} {rmsd_sd:.2f} {rmsdemin_mean:.2f} {rmsdemin_sd:.2f} {size} {total_mean:.2f} {total_sd:.2f} {vdw_mean:.2f} {vdw_sd:.2f} {elec_mean:.2f} {elec_sd:.2f} {bsa_mean:.2f} {bsa_sd:.2f}\n')
			f.close()

		return output_f

	def calculate_contacts(self):
		""" Calculate atomic contacts """
		# SORT BY HADDOCK-SCORE
		pdb_list = [(e, self.structure_dic[e]['haddock-score']) for e in self.structure_dic]
		sorted_pdb_list = sorted(pdb_list, key=lambda x: x[1])

		with open('file.list', 'w') as fl_fh:
			for p in sorted_pdb_list:
				fl_fh.write(f'{p[0]}\t{{ {p[1]:.2f} }}\n')
		fl_fh.close()

		file_nam = [p[0] for p in sorted_pdb_list]

		# with open('file.nam', 'w') as fh:
		# 	for p in file_nam:
		# 		fh.write(f'{p}\n')
		# fh.close()

		con = Contacts()
		con.calculate_contacts(file_nam, cutoff=5.0)

		self.con_list = con.contact_file_list

	def calc_fcc_matrix(self, ignore_chain=False):
		""" Calculate the FCC matrix (extracted and adapted from calc_fcc_matrix.py """
		fcc_matrix_outf = 'fcc.matrix'

		if not self.con_list:
			self.calculate_contacts()

		# contacts = []
		# for con_f in self.con_list:
		# 	t = []
		# 	with open(con_f) as f:
		# 		for l in f.readlines():
		# 			t.append(int(l))
		# 	f.close()
		# 	contacts.append(set(t))

		if ignore_chain:
			contacts = [[int(l[0:5] + l[6:-1]) for l in open(f)] for f in self.con_list if f.strip()]
		else:
			contacts = [set([int(l) for l in open(f)]) for f in self.con_list if f.strip()]

		with open(fcc_matrix_outf, 'w') as fcc_fh:
			# get the pairwise matrix
			for i in range(len(contacts)):
				for j in range(i + 1, len(contacts)):
					con_a = contacts[i]
					con_b = contacts[j]

					cc = float(len(con_a.intersection(con_b)))
					# cc_v = float(len(con_b.intersection(con_a)))
					# def calculate_fcc(listA, listB):
					# 	"""
					# 	Calculates the fraction of common elements between two lists
					# 	taking into account chain IDs
					# 	"""
					#
					# 	cc = len(listA.intersection(listB))
					# 	cc_v = len(listB.intersection(listA))
					#
					# 	return (cc, cc_v)

					try:
						fcc, fcc_v = cc * 1.0 / len(con_a), cc * 1.0 / len(con_b)
					except ZeroDivisionError:
						# Either con_a or con_b == 0, which means that there are no contacts between molecules.
						#  This should not happen in docking, but is relevant for scoring
						fcc = .0
						fcc_v = .0

					self.fcc_matrix.append(f'{i + 1} {j + 1} {fcc:.3f} {fcc_v:.3f}')

					fcc_fh.write(f'{i + 1} {j + 1} {fcc:.3f} {fcc_v:.3f}\n')

		fcc_fh.close()
		return fcc_matrix_outf

	def fetch_lowest(self):
		score_list = []
		for pdb in self.structure_dic:
			score_list.append((pdb, self.structure_dic[pdb]['haddock-score']))
		sorted_score_list = sorted(score_list, key=lambda x: x[1])
		lowest_pdb = sorted_score_list[0][0]
		lowest_score = sorted_score_list[0][1]
		return lowest_pdb, lowest_score

	def run_dockq(self, ref):

		dockq_exec = ini.get('third_party', 'dockq_exe')

		if shutil.which(dockq_exec):
			if ref == 'None':
				print(f'++ Cannot run DockQ, no reference has been set')
				return

			print('\n+ Running DockQ')

			if ref == 'lowest':
				reference_pdb, reference_score = self.fetch_lowest()
				total = len(self.structure_dic)
				print(f'++ Using lowest haddock score as reference {reference_pdb} = {reference_score:.3f}, n = {total}')

			else:
				reference_pdb = ref

			for i, pdb in enumerate(self.structure_dic):
				result_dic = dockq(reference_pdb, pdb, dockq_exec)
				for k in result_dic:
					self.structure_dic[pdb][k] = result_dic[k]
		else:
			print('\n+ DockQ not correctly configured in haddock3.ini')
			return

	def run_fastcontact(self):
		""" Run fastcontact on all scored PDBs """

		fastcontact_exe = ini.get('third_party', 'fastcontact_exe')

		if shutil.which(fastcontact_exe):
			print('\n+ Running FASTCONTACT')

			for pdb in self.structure_dic:
				fast_elec, fast_desol = fastcontact(pdb, fastcontact_exe)
				self.structure_dic[pdb]['fastelec'] = fast_elec
				self.structure_dic[pdb]['fastdesol'] = fast_desol
		else:
			print('\n+ FASTCONTACT not correctly configured in haddock3.ini')
			return

	def run_dfire(self):
		""" Run dfire on all scored PDBs """

		dcomplex_exe = ini.get('third_party', 'dcomplex_exe')
		dcomplex_lib = ini.get('third_party', 'dcomplex_lib')
		dcomplex_dat = ini.get('third_party', 'dcomplex_dat')

		if shutil.which(dcomplex_exe) and os.path.exists(dcomplex_lib) and os.path.exists(dcomplex_dat):
			print('\n+ Running DFIRE')
			current_path = os.getcwd()
			shutil.copy(dcomplex_lib, current_path)
			shutil.copy(dcomplex_dat, current_path)

			for pdb in self.structure_dic:
				d_binding, d_score = dfire(pdb, dcomplex_exe)
				self.structure_dic[pdb]['dfire-ebinding'] = d_binding
				self.structure_dic[pdb]['dfire-score'] = d_score

			lib_name = dcomplex_lib.split('/')[-1]
			dat_name = dcomplex_dat.split('/')[-1]
			os.remove(lib_name)
			os.remove(dat_name)
		else:
			print('\n+ DFIRE not correctly configured in haddock3.ini')
			return

	def output(self):

		# Single Structure
		print('\n+ Saving single-structure analysis')

		# sort by haddock score!
		score_list = [(pdb, self.structure_dic[pdb]['haddock-score']) for pdb in self.structure_dic]
		sorted_score_list = sorted(score_list, key=lambda x: x[1])
		k = sorted_score_list[0][0]
		header = 'model ranking ' + ' '.join(list(self.structure_dic[k])) + '\n'

		with open(f'analysis/{self.ss_output_f}.stats', 'w') as out:
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

		# Cluster based - improve ASAP #======================================================#
		print('\n+ Saving cluster analysis')

		cluster_params = list(set([c.split('cluster-')[-1].split('_')[0] for c in self.structure_dic[k] if 'cluster' in c]))
		cluster_dic = {}
		for param in cluster_params:
			cluster_dic[param] = []
			for pdb in self.structure_dic:

				cluster_id = self.structure_dic[pdb][f'cluster-{param}_name']
				cluster_ranking = self.structure_dic[pdb][f'cluster-{param}_overall_ranking']
				cluster_internal_ranking = self.structure_dic[pdb][f'cluster-{param}_internal_ranking']

				if not math.isnan(cluster_id):

					cluster_dic[param].append( (cluster_ranking, cluster_id, cluster_internal_ranking, pdb, self.structure_dic[pdb]) )

		# Sort cluster_dic according to cluster_internal_ranking
		header = 'cluster size haddock-score sd bsa sd Edesolv sd Evdw sd fc-elec sd fc-desol sd Dfire-Ebinding sd dfire-score sd\n'
		for param in cluster_dic:

			cluster_ranking_list = list(set([ e[0] for e in cluster_dic[param] ]))
			cluster_ranking_list.sort()

			tbw = ''
			tbw_4 = ''
			tbw_2 = ''

			for rank in cluster_ranking_list:

				# get all elements of this rank
				elements = [e for e in cluster_dic[param] if e[0] == rank]
				cluster_id = elements[0][1]

				element_list = [(j[2], j[3], j[4]) for j in elements]
				sorted_element_list = sorted(element_list, key=lambda x: x[0])

				# was the cluster center added to cluster_dic? if so remove!
				# del sorted_element_list[0]

				# start!

				# Haddock Score
				haddock_score_l = [j[2]['haddock-score'] for j in sorted_element_list]

				haddock_score_mean = statistics.mean(haddock_score_l)
				haddock_score_sd = statistics.stdev(haddock_score_l)

				haddock_score_mean_4 = statistics.mean(haddock_score_l[:4])
				haddock_score_sd_4 = statistics.stdev(haddock_score_l[:4])

				haddock_score_mean_2 = statistics.mean(haddock_score_l[:2])
				haddock_score_sd_2 = statistics.stdev(haddock_score_l[:2])

				# BSA
				bsa_l = [j[2]['bsa'] for j in sorted_element_list]
				bsa_mean = statistics.mean(bsa_l)
				bsa_sd = statistics.stdev(bsa_l)

				bsa_mean_4 = statistics.mean(bsa_l[:4])
				bsa_sd_4 = statistics.stdev(bsa_l[:4])

				bsa_mean_2 = statistics.mean(bsa_l[:2])
				bsa_sd_2 = statistics.stdev(bsa_l[:2])

				# Desolv
				desolv_l = [j[2]['desolv'] for j in sorted_element_list]
				desolv_mean = statistics.mean(desolv_l)
				desolv_sd = statistics.stdev(desolv_l)

				desolv_mean_4 = statistics.mean(desolv_l[:4])
				desolv_sd_4 = statistics.stdev(desolv_l[:4])

				desolv_mean_2 = statistics.mean(desolv_l[:2])
				desolv_sd_2 = statistics.stdev(desolv_l[:2])

				# VdW
				vdw_l = [j[2]['vdw'] for j in sorted_element_list]
				vdw_mean = statistics.mean(vdw_l)
				vdw_sd = statistics.stdev(vdw_l)

				vdw_mean_4 = statistics.mean(vdw_l[:4])
				vdw_sd_4 = statistics.stdev(vdw_l[:4])

				vdw_mean_2 = statistics.mean(vdw_l[:2])
				vdw_sd_2 = statistics.stdev(vdw_l[:2])

				# Size
				size = len(haddock_score_l)

				try:
					fastelec_l = [j[2]['fastelec'] for j in sorted_element_list]
					fastelec_mean = statistics.mean(fastelec_l)
					fastelec_sd = statistics.stdev(fastelec_l)

					fastelec_mean_4 = statistics.mean(fastelec_l[:4])
					fastelec_sd_4 = statistics.stdev(fastelec_l[:4])

					fastelec_mean_2 = statistics.mean(fastelec_l[:2])
					fastelec_sd_2 = statistics.stdev(fastelec_l[:2])

				except:
					fastelec_mean = .0
					fastelec_sd = .0

					fastelec_mean_4 = .0
					fastelec_sd_4 = .0

					fastelec_mean_2 = .0
					fastelec_sd_2 = .0

				try:
					fastdesol_l = [j[2]['fastdesol'] for j in sorted_element_list]
					fastdesol_mean = statistics.mean(fastdesol_l)
					fastdesol_sd = statistics.stdev(fastdesol_l)

					fastdesol_mean_4 = statistics.mean(fastdesol_l[:4])
					fastdesol_sd_4 = statistics.stdev(fastdesol_l[:4])

					fastdesol_mean_2 = statistics.mean(fastdesol_l[:2])
					fastdesol_sd_2 = statistics.stdev(fastdesol_l[:2])

				except:
					fastdesol_mean = .0
					fastdesol_sd = .0

					fastdesol_mean_4 = .0
					fastdesol_sd_4 = .0

					fastdesol_mean_2 = .0
					fastdesol_sd_2 = .0

				try:
					dfire_ebinding_l = [j[2]['dfire-ebinding'] for j in sorted_element_list]
					dfire_ebinding_mean = statistics.mean(dfire_ebinding_l)
					dfire_ebinding_sd = statistics.stdev(dfire_ebinding_l)

					dfire_ebinding_mean_4 = statistics.mean(dfire_ebinding_l[:4])
					dfire_ebinding_sd_4 = statistics.stdev(dfire_ebinding_l[:4])

					dfire_ebinding_mean_2 = statistics.mean(dfire_ebinding_l[:2])
					dfire_ebinding_sd_2 = statistics.stdev(dfire_ebinding_l[:2])

					dfire_score_l = [j[2]['dfire-score'] for j in sorted_element_list]
					dfire_score_mean = statistics.mean(dfire_score_l)
					dfire_score_sd = statistics.stdev(dfire_score_l)

					dfire_score_mean_4 = statistics.mean(dfire_score_l[:4])
					dfire_score_sd_4 = statistics.stdev(dfire_score_l[:4])

					dfire_score_mean_2 = statistics.mean(dfire_score_l[:2])
					dfire_score_sd_2 = statistics.stdev(dfire_score_l[:2])

				except:
					dfire_ebinding_mean = .0
					dfire_ebinding_sd = .0

					dfire_ebinding_mean_4 = .0
					dfire_ebinding_sd_4 = .0

					dfire_ebinding_mean_2 = .0
					dfire_ebinding_sd_2 = .0

					dfire_score_mean = .0
					dfire_score_sd = .0

					dfire_score_mean_4 = .0
					dfire_score_sd_4 = .0

					dfire_score_mean_2 = .0
					dfire_score_sd_2 = .0

				tbw += f'clust{cluster_id} {size} {haddock_score_mean:.2f} {haddock_score_sd:.2f} {bsa_mean:.2f} {bsa_sd:.2f} {desolv_mean:.2f} {desolv_sd:.2f} {vdw_mean:.2f} {vdw_sd:.2f} {fastelec_mean:.2f} {fastelec_sd:.2f} {fastdesol_mean:.2f} {fastdesol_sd:.2f} {dfire_ebinding_mean:.2f} {dfire_ebinding_sd:.2f} {dfire_score_mean:.2f} {dfire_score_sd:.2f}\n'
				tbw_4 += f'clust{cluster_id} {size} {haddock_score_mean_4:.2f} {haddock_score_sd_4:.2f} {bsa_mean_4:.2f} {bsa_sd_4:.2f} {desolv_mean_4:.2f} {desolv_sd_4:.2f} {vdw_mean_4:.2f} {vdw_sd_4:.2f} {fastelec_mean_4:.2f} {fastelec_sd_4:.2f} {fastdesol_mean_4:.2f} {fastdesol_sd_4:.2f} {dfire_ebinding_mean_4:.2f} {dfire_ebinding_sd_4:.2f} {dfire_score_mean_4:.2f} {dfire_score_sd_4:.2f}\n'
				tbw_2 += f'clust{cluster_id} {size} {haddock_score_mean_2:.2f} {haddock_score_sd_2:.2f} {bsa_mean_2:.2f} {bsa_sd_2:.2f} {desolv_mean_2:.2f} {desolv_sd_2:.2f} {vdw_mean_2:.2f} {vdw_sd_2:.2f} {fastelec_mean_2:.2f} {fastelec_sd_2:.2f} {fastdesol_mean_2:.2f} {fastdesol_sd_2:.2f} {dfire_ebinding_mean_2:.2f} {dfire_ebinding_sd_2:.2f} {dfire_score_mean_2:.2f} {dfire_score_sd_2:.2f}\n'

			# with open(f'{self.cluster_output_f}_{param}.stats', 'w') as f:
			# 	f.write(header)
			# 	f.write(tbw)
			# f.close()

			# with open(f'{self.cluster_output_f}_{param}_best4.stats', 'w') as f:
			# 	f.write(header)
			# 	f.write(tbw_4)
			# f.close()

			# with open(f'{self.cluster_output_f}_{param}_best2.stats', 'w') as f:
			# 	f.write(header)
			# 	f.write(tbw_2)
			# f.close()

	def match_renumber(self, reference_pdb):
		""" Match the chains and renumber the structures according to a reference PDB """

		clustalo_exe = ini.get('third_party', 'clustalo_exe')

		if not shutil.which(clustalo_exe):
			print('\n+ ClustalO not correctly configured in haddock3.ini')
			print('+ WARNING: matching not possible!')
			return False
		else:
			print('\n+ Running automated chain matching and renumbering')
			print('+ WARNING: Use with caution, some residues could be deleted')

		if reference_pdb == 'lowest':
			reference_pdb, reference_score = self.fetch_lowest()

		if ' ' not in PDB.identify_chains(reference_pdb):
			reference_pdb = PDB.fix_id(reference_pdb, priority='chain')
		if ' ' not in PDB.identify_chainseg(reference_pdb):
			reference_pdb = PDB.fix_id(reference_pdb, priority='seg')

		reference_seq_dic = PDB.load_seq(reference_pdb)
		reference_chains = PDB.identify_chains(reference_pdb)
		reference_chains.sort()

		pdb_list = list(self.structure_dic.keys())
		pdb_list.sort()

		for pdb in pdb_list:
			# match the chains with sequence alignment
			target_seq_dic = PDB.load_seq(pdb)
			pdb = PDB.fix_id(pdb, priority='seg')

			# Get what chains are present in the target
			target_chains = PDB.identify_chains(pdb)
			target_chains.sort()

			# Do a combinatorial alignment to check which chains match better
			identity_dic = {}
			for ref_chain, target_chain in itertools.product(reference_chains, target_chains):
				ref_seq = ''.join(list(reference_seq_dic[ref_chain].values()))
				target_seq = ''.join(list(target_seq_dic[target_chain].values()))
				open('seq.fasta', 'w').write(f'>ref\n{ref_seq}\n>target\n{target_seq}\n')
				cmd = f'{clustalo_exe} -i seq.fasta --outfmt=clu --resno --wrap=9000 --force'
				p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out = p.communicate()
				os.remove('seq.fasta')
				aln_data = out[0].decode('utf-8').split()
				ref_aln = aln_data[6]
				target_aln = aln_data[9]
				counter_a = 0
				counter_b = 0
				numbering_dic = {}
				for i in range(len(ref_aln)):
					ref_char = ref_aln[i]
					target_char = target_aln[i]
					ref_resnum = list(reference_seq_dic[ref_chain])[counter_a]
					try:
						target_resnum = list(target_seq_dic[target_chain])[counter_b]
					except IndexError:
						# Target sequence exhausted, ignore
						target_resnum = '-'
					# print(ref_char, ref_resnum, target_char, target_resnum)
					if '-' not in ref_char:
						counter_a += 1
					if '-' not in target_char:
						counter_b += 1
					if '-' not in ref_char and '-' not in target_char:
						numbering_dic[target_resnum] = ref_resnum
				identity = out[0].decode('utf-8').count('*') / float(len(ref_seq))
				coverage = len(numbering_dic) / len(ref_aln)
				# print(ref_chain, target_chain, identity, coverage)
				# print(f'>R:{ref_chain}\n{ref_aln}')
				# print(f'>T:{target_chain}\n{target_aln}')
				try:
					identity_dic[ref_chain].append((target_chain, identity, coverage, numbering_dic))
				except KeyError:
					identity_dic[ref_chain] = [(target_chain, identity, coverage,  numbering_dic)]
			# print('#########\n')

			# do the renumbering
			for i, ref_c in enumerate(reference_chains):
				target_info_list = [(v[0], v[1], v[2]) for v in identity_dic[ref_c]]
				# sort by identity and coverage
				sorted_target_list = sorted(target_info_list, key=lambda x: (-x[2], x[1]))
				# create a catalog with possible numbering references
				numbering_dic_catalog = dict([(v[0], v[3]) for v in identity_dic[ref_c]])

				if len(set([ e[1] for e in sorted_target_list])) == 1:
					# this is a homo-something, match is sequentialy
					selected_chain = target_chains[i]
				else:
					# get the highest identity/coverage
					selected_chain = sorted_target_list[0][0]

				# select the correct numbering dictionary
				selected_numbering_dic = numbering_dic_catalog[selected_chain]
				# just for readability:
				old_chain = selected_chain
				new_chain = ref_c

				# replace the target chain (old) with the same observed in the reference (new)
				chain_matched_pdb = PDB.replace_chain(pdb, old_chain, new_chain)

				# renumber!
				#  Note, if residue is present in target and not in reference it will be DELETED, use with caution
				renumbered_pdb = PDB.renumber(chain_matched_pdb, selected_numbering_dic, new_chain, overwrite=True)

		return True

	# =================================================================================================================#


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
