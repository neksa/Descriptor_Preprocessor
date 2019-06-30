import os
import config
import pickle

from collections import defaultdict, Counter
import math
from operator import itemgetter
import os
import pickle as pkl

import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import numpy as npno

# import seaborn
import config
from utils import seq_logo
import numpy as np

def plot_dihedral_for_diff_res(descr_full, sno_position):
    descr_full = descr_full.groupby('relative_sno')
    descr_for_sno = descr_full.get_group(sno_position)
    descr_for_sno_grouped = descr_for_sno.groupby('res')
    num_plots = len(descr_for_sno_grouped)
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)
    for ax_count, (res, descr_for_res) \
            in enumerate(descr_for_sno_grouped, start=1):
        phi = descr_for_res.phi.values
        psi = descr_for_res.psi.values
        ax = plt.subplot(num_row, num_col, ax_count)
        ax.set_title(str(res))
        ax.scatter(phi, psi, marker='.', s=1)
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])

#pylint: disable=invalid-name
def plot_CA(descr_full):
    descr_full_grouped = descr_full.groupby('filename')
    for filename, df_per_file in descr_full_grouped:
        plt.figure()
        plt.suptitle(filename)
        df_per_file = df_per_file.groupby(['cid', 'seq_marker'])
        num_plots = len(df_per_file)
        num_col = math.ceil(math.sqrt(num_plots))
        num_row = math.ceil(num_plots / num_col)
        for ax_count, ((cid, seq_marker), df_per_set) in \
                enumerate(df_per_file, start=1):
            # subplot indices start from 1
            ax = plt.subplot(num_row, num_col, ax_count, projection='3d')
            ax.set_title(str(seq_marker) + " : " + str(cid))

            snos, CA_coords = df_per_set.sno, df_per_set.CA
            for i in CA_coords:
                assert len(i) == 3
            xs, ys, zs = list(zip(*CA_coords))

            # c for color distribution
            ax.scatter3D(xs, ys, zs, c=range(len(xs)), marker='o',
                         cmap='autumn')
            ax.plot3D(xs, ys, zs, 'k', linewidth=0.5)

            # Mark out 0th position on plot
            for i, sno in enumerate(snos):
                if sno == seq_marker:
                    ax.text(xs[i], ys[i], zs[i], sno)

def _for_sort(x):
    return -x[1], x[0]

def plot_signature_logo(descr_full, num_res_considered=4):
    descr_filtered = descr_full.filter(['relative_sno', 'res'])
    min_sno = min(descr_filtered.relative_sno)
    descr_sno_grouped = descr_filtered.groupby('relative_sno')
    to_logo = [[] for __ in range(len(descr_sno_grouped))]
    for i, (_, descr_per_sno) in enumerate(descr_sno_grouped):
        # key is sno
        res_unsorted = descr_per_sno.res.values
        res_unique = list(np.unique(res_unsorted, return_counts=True))
        num_seqs = sum(res_unique[1])
        res_sorted = list(zip(*sorted(zip(*res_unique), key=_for_sort)))
        res_names, res_counts = res_sorted[0][:num_res_considered], \
                                res_sorted[1][:num_res_considered]
        res_names = res_names[::-1]
        res_counts = res_counts[::-1]
        res_percents = list([i/num_seqs for i in res_counts])
        for name, percent in zip(res_names, res_percents):
            to_logo[i].append((name, percent))
    seq_logo.Logo(to_logo, min_sno)

ref_ = os.path.join(config.ROOT, 'tests', 'data', 'ref')
prosite_output = os.path.join(ref_, 'motif_pos_prosite.pkl')

with open(prosite_output, 'rb') as file:
    motif_pos = pickle.load(file)
from descr import descr_main, loaders

pdb_info_data_map = loaders.load_pdb_info(motif_pos)
descrs = descr_main.calculate(pdb_info_data_map)
print(descrs)
# with open("tmp.pkl", 'rb') as file:
#     descrs = pickle.load(file)
print(descrs.columns)
print(set(descrs.filename))
print(set(descrs.relative_sno))
print(len(descrs))

plot_signature_logo(descrs)
plt.show()