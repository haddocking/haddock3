# compare scores (:
import matplotlib.pyplot as plt
import pandas as pd

ens_f = 'all-runs-combined.pdb'
file_l = 'run-all-runs-combined/file.list'
golden_f = '/home/abonvin/target161-analysis/all-runs-combined/file.list'

ens_d = {}
file_dic = {}
with open(file_l, 'r') as file_f:
	for l in file_f.readlines():
		data = l.split()
		model_id = int(data[0].split('_')[-1].split('.pdb')[0])
		score = float(data[-2])
		file_dic[model_id] = score


with open(ens_f, 'r') as ens_f:
	for l in ens_f.readlines():
		if l.startswith('REMARK'):
			data = l.split()
			model_id = int(data[2])
			full_name = data[4]
			score = file_dic[model_id]
			if score <= 200:
				ens_d[full_name] = [score]

with open(golden_f, 'r') as g_f:
	for l in g_f.readlines():
		data = l.split()
		golden_name = ''.join(data[0].split('_conv'))
		score = float(data[-2])
		try:
			ens_d[golden_name].append(score)
		except:
			pass

df = pd.DataFrame.from_dict(ens_d, orient='index', columns=['v3.0', 'v2.4'])
df.plot.scatter(x='v3.0', y='v2.4', title='Single structure')
plt.savefig('ss-scores.png')

# compare top-cluster scores?

cluster_f = '/home/rodrigo/projects/scoring/T161/rescoring/run-all-runs-combined/analysis/fcc_0.6-4.0.stats_best4'
# load and sort by haddockscore

cluster_scores = []
with open(cluster_f, 'r') as fh:
	for l in fh.readlines():
		if 'f' in l[0]:
			score = float(l.split()[1])
			cluster_scores.append(score)

cluster_scores.sort()


# already sorted
golden_f = '/home/abonvin/target161-analysis/all-runs-combined/new-clusters-min4/statistics'
golden_scores = []
with open(golden_f, 'r') as fh:
	for l in fh.readlines():
		if 'c' in l[0]:
			score = float(l.split()[4])
			golden_scores.append(score)

if len(cluster_scores) > len(golden_scores):
	cutoff = len(golden_scores)
else:
	cutoff = len(cluster_scores)

df = pd.DataFrame.from_records(list(zip(cluster_scores[:cutoff], golden_scores[:cutoff])), columns=['v3.0','v2.4'])
df.plot.scatter(x='v3.0', y='v2.4', title='Top clusters average')
plt.savefig('clt-scores.png')

