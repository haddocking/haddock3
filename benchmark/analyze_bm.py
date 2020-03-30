# This should be ported to a notebook
import re
import math
import argparse
import glob
import pandas as pd
import numpy as np
# parse the BM5 run and prepare an analysis friendly file

# expected data format is:
#       PDB1      PDB2       PDB3    ...  PDB50     PDB51     PDB52
# 0    2.926811  0.326559  1.833239  ...  1.540668  3.324400  2.635251
# 1    1.258266  1.181834  2.798129  ...  3.501930  2.065754  1.658211
# 2    3.813568  2.633314  0.049096  ...  0.565649  2.561827  0.389702


def parse_ss_data(ss_f):
    data_l = []
    name = re.findall(r'(\w*)(\/|.)run', ss_f)[0]
    with open(ss_f) as fh:
        data = fh.readlines()
        header = data[0].split()
        irms_index = [i for i, j in enumerate(header) if 'irms' in j or 'i-RMSD' in j]
        if len(irms_index) > 1:
            print('+ ERROR, more than one interface detected, this is not yet supported!')
            exit()
        else:
            irms_index = irms_index[0]

        for entry in data[1:]:
            irmsd = float(entry.split()[irms_index])
            data_l.append(irmsd)

    return name, data_l


def parse_clt_data(ss_f):
    data_l = []
    name = re.findall(r'(\w*)(\/|.)run', ss_f)[0]
    with open(ss_f) as fh:
        clt_dic = {}
        data = fh.readlines()
        header = data[0].split()
        # clt_name_index = [i for i, j in enumerate(header) if 'cluster-0.6-4_name' in j][0]
        clt_ranking_index = [i for i, j in enumerate(header) if 'cluster-0.6-4_overall_ranking' in j][0]
        irms_index = [i for i, j in enumerate(header) if 'irms' in j or 'i-RMSD' in j]
        if len(irms_index) > 1:
            print('+ ERROR, more than one interface detected, this is not yet supported!')
            exit()
        else:
            irms_index = irms_index[0]
        for entry in data[1:]:
            # clt_name = int(entry.split()[clt_name_index])
            clt_ranking = float(entry.split()[clt_ranking_index])
            irmsd = float(entry.split()[irms_index])
            if not math.isnan(clt_ranking):
                if clt_ranking not in clt_dic:
                    clt_dic[clt_ranking] = []
                clt_dic[clt_ranking].append(irmsd)

        for c in clt_dic:
            min_irmsd = min(clt_dic[c])
            data_l.append(min_irmsd)

    return name, data_l


def analyze(ss_l, ss=None, clt=None):

    bm_dic = {}
    top_range = []
    identifier = None
    values = None
    for i, stat_f in enumerate(ss_l):
        if ss:
            identifier, values = parse_ss_data(stat_f)
            top_range = [1, 5, 10, 50, 100, 200]
        if clt:
            identifier, values = parse_clt_data(stat_f)
            top_range = [1, 2, 5, 10]

        bm_dic[identifier] = values

    # create as index orientation to account for arrays with different lenghts (aka clusters)
    #  its slow but gets the job done (:
    i_df = pd.DataFrame.from_dict(bm_dic, orient='index')
    df = i_df.transpose()

    total = len(df.columns)
    bm_dic = get_stars(df, top_range)

    bm_df = pd.DataFrame.from_dict(bm_dic)
    per_bm_df = ((bm_df / total) * 100).round(2)

    bm_md = bm_df.to_markdown().split('\n')
    bm_md[0] = '|   ' + ''.join([f'| Top{i} ' for i in top_range]) + '|'
    bm_md = '\n'.join(bm_md)
    bm_md = re.sub(r'\|(\s\s0\s)', '| #success *', bm_md)
    bm_md = re.sub(r'\|(\s\s1\s)', '| #success **', bm_md)
    bm_md = re.sub(r'\|(\s\s2\s)', '| #success ***', bm_md)

    per_bm_md = per_bm_df.to_markdown().split('\n')
    per_bm_md[0] = '|   ' + ''.join([f'| Top{i} ' for i in top_range]) + '|'
    per_bm_md = '\n'.join(per_bm_md)
    per_bm_md = re.sub(r'\|(\s\s0\s)', '| %success *', per_bm_md)
    per_bm_md = re.sub(r'\|(\s\s1\s)', '| %success **', per_bm_md)
    per_bm_md = re.sub(r'\|(\s\s2\s)', '| %success ***', per_bm_md)

    return f'\nTOTAL: {total}\n{bm_md}\n\n{per_bm_md}'


def get_stars(data_frame, top_r):
    d = {}
    for top in top_r:
        subset_df = data_frame[:top].min()
        high = (subset_df <= 1).sum()
        # medium = ((subset_df > 1) & (subset_df <= 2)).sum()
        medium = (subset_df <= 2).sum() - 1
        # low = ((subset_df > 2) & (subset_df <= 4)).sum()
        low = (subset_df <= 4).sum() - 1
        d[top] = low, medium, high
    return d


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare SS data')
    parser.add_argument("benchmark_path",
                        help="Location of BM5")

    parser.add_argument("--haddock3",
                        action='store_true',
                        default=False)

    parser.add_argument("--haddock24",
                        action='store_true',
                        default=False)

    args = parser.parse_args()

    path = args.benchmark_path

    ss_list = []
    stats_path = ''
    if args.haddock3:
        stats_path = f'{path}/*/run-ti/ss.stats'
    elif args.haddock24:
        stats_path = f'{path}/*/*run1-ti_structures-water.stats'

    ss_list = glob.glob(stats_path)
    ss_ana = analyze(ss_list, ss=True)

    cluster_ana = None
    if args.haddock3:
        cluster_ana = analyze(ss_list, clt=True)
    elif args.haddock24:
        cluster_ana = 'Not supported'

    print('\n++ Single-structure', ss_ana)
    print('\n++ Cluster', cluster_ana)

