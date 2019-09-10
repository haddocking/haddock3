import glob


def read_table(data_f):
	d = {}
	with open(data_f) as f:
		for rank, line in enumerate(f.readlines()[1:]):
			pdb, score, irmsd, lrmsd, fnat = line.split('\t')
			d[rank+1] = pdb, float(score), float(irmsd), float(lrmsd), float(fnat)
	return d


def calc_capri_stats(data_dic, top):
	high = 0
	medium = 0
	acceptable = 0
	for i, rank in enumerate(data_dic):
		if i == top:
			break
		pdb, score, irmsd, lrmsd, fnat = data_dic[rank]

		if (fnat >= 0.5) and (lrmsd <= 1.0) or (irmsd <= 1.0):
			high += 1
		elif (fnat >= 0.3) and (lrmsd <= 5.0) or (irmsd <= 2.0):
			medium += 1
		elif (fnat >= 0.1) and (lrmsd <= 10.0) or (irmsd <= 4.0):
			acceptable += 1

	return high, medium, acceptable

# tbw = 'target\thigh\tmedium\tacceptable\thigh\tmedium\tacceptable\n'
# tbw = 'target\tmedium\tacceptable\tmedium\tacceptable\n'
for folder in glob.glob('*'):

	haddock24_data = read_table(f'{folder}/it0.haddock24.dat')
	haddock3_data = read_table(f'{folder}/it0.haddock3.dat')

	# number of acceptable / medium / high

	haddock24_stats = calc_capri_stats(haddock24_data, 400)
	haddock3_stats = calc_capri_stats(haddock3_data, 400)

	print(f'{folder}\t{haddock24_stats[0]}/*\t\t{haddock24_stats[1]}/**\t{haddock24_stats[2]}/***')
	print(f'\t\t{haddock3_stats[0]}/*\t\t{haddock3_stats[1]}/**\t{haddock3_stats[2]}/***')
	print('\n')


# print(tbw)