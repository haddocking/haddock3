# read all benchmark data and generate some plots
import argparse
import collections
import re
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_haddockscore(ss_d, irmsd_load=False):
    c = 0
    d = {}
    i_d = {}
    for k in ss_d:
        stat_f = ss_d[k]
        # name = re.findall(r'(\w*)(\/|.)run', stat_f)[0]
        with open(stat_f) as fh:
            data = fh.readlines()
            header = data[0].split()
            hs_index = [i for i, j in enumerate(header) if 'haddock-score' in j][0]
            if irmsd_load:
                irms_index = [i for i, j in enumerate(header) if 'irms' in j or 'i-RMSD' in j][0]
            for entry in data[1:]:
                hs = float(entry.split()[hs_index])
                if irmsd_load:
                    irmsd =float(entry.split()[irms_index])
                    i_d[c] = irmsd
                d[c] = hs
                c += 1
    if irmsd_load:
        return d, i_d
    else:
        return d


def get_names(ss_l):
    d = {}
    for ss in ss_l:
        name = re.findall(r'\/(\w*)(\/|).run', ss)[-1][0]
        d[name] = ss
    return d



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare SS data')

    parser.add_argument("bm24_path")
    parser.add_argument("bm30_path")
    args = parser.parse_args()

    bm24_stats = f'{args.bm24_path}/*/*run1-ti_structures-water.stats'
    bm30_stats = f'{args.bm30_path}/*/run-ti/ss.stats'

    # make sure its paired!
    bm24_d = get_names(glob.glob(bm24_stats))
    bm30_d = get_names(glob.glob(bm30_stats))

    bm24_od = collections.OrderedDict(sorted(bm24_d.items()))
    bm30_od = collections.OrderedDict(sorted(bm30_d.items()))

    # [e for e in bm24_od] == [f for f in bm30_od]

    bm24_hs_dic, bm24_irmsd_dic = load_haddockscore(bm24_od, irmsd_load=True)
    bm30_hs_dic, bm30_irmsd_dic = load_haddockscore(bm30_od, irmsd_load=True)

    df = pd.DataFrame([bm24_hs_dic, bm24_irmsd_dic, bm30_hs_dic, bm30_irmsd_dic])
    t_df = df.transpose()
    t_df = t_df.rename(columns={0:'v2.4', 1:'v2.4-irmsd', 2: 'v3.0', 3: 'v3.0-irmsd'})

    # plots #####

    # kdeplot
    f, axes = plt.subplots(1, 2,figsize=(12,6))
    f.tight_layout(pad=3.0)
    p = sns.kdeplot(t_df['v2.4'], t_df['v2.4-irmsd'], cmap="Blues", shade=True, ax=axes[0], n_levels=15)
    p.set(xlim=(-300, 100))
    p.set(ylim=(0, 17.5))
    p.set(xlabel='Haddock-score', ylabel='i-RMSD',title='v2.4')
    p = sns.kdeplot(t_df['v3.0'], t_df['v3.0-irmsd'], cmap="Oranges", shade=True, ax=axes[1], n_levels=15)
    p.set(xlim=(-300, 100))
    p.set(ylim=(0, 17.5))
    p.set(xlabel='Haddock-score', ylabel='i-RMSD',title='v3.0')
    plt.savefig('kde.png')
    plt.clf()

    # jointplot
    p = sns.jointplot(x="v2.4", y="v3.0", data=t_df, kind="kde",xlim=(-300, 100), ylim=(-300, 100), color="Blue")
    plt.savefig('jointplot.png')
    plt.clf()

    # histogram
    sns.distplot(t_df['v2.4'], kde=True, norm_hist=True, hist=True, bins=100)
    p = sns.distplot(t_df['v3.0'], kde=True, norm_hist=True, hist=True,bins=100)
    p.set(xlim=(-300, 500))
    plt.savefig('histogram.png')
    plt.clf()
