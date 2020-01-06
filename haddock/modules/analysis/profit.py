# from haddock.workflows.scoring.config import load_parameters
#
# param_dic = load_parameters()
# profit_exe = param_dic['third-party']['profit_exe']
#
#
# def profit_rmsd(ordered_dic):
# 	# profit_exe = '/programs/i386-mac/profit/3.1/profit'
# 	rmsd_d = {}
#
# 	# sort by haddock score
# 	hs_list = [(k, struc_dic[k]['haddock-score']) for k in struc_dic]
# 	sorted_hs_list = sorted(hs_list, key=lambda x: (-x[1], x[0]))
# 	sorted_hs_list.reverse()
# 	sorted_pdb_list = [e[0] for e in sorted_hs_list]
#
# 	refe = sorted_pdb_list[0]
# 	for mobi in sorted_pdb_list:
# 		cmd = 'refe {}\nmobi {}\nignore missing\natoms CA\nzone A*\nfit\nrzone B*\nquit'.format(refe, mobi)
# 		output = os.popen('echo "{}" | {}'.format(cmd, profit_exe))
# 		result = [l for l in output if 'RMS:' in l][-1]
# 		rmsd = float(result.split()[-1])
# 		rmsd_d[mobi] = rmsd
#
# 	return rmsd_d