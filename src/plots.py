from collections import defaultdict, Counter
import math
from operator import itemgetter
import os
import pickle as pkl

import matplotlib.pyplot as plt
import numpy as np

from config import paths
from utils import seq_logo, generic

def plot_dihedral_for_diff_res(descr_full, sno_position):
    descr_full = descr_full.groupby('relative_sno')
    descr_for_sno = descr_full.get_group(sno_position)
    descr_for_sno_grouped = descr_for_sno.groupby('res')
    num_plots = len(descr_for_sno_grouped)
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)
    for ax_count, (res, descr_for_res) in enumerate(descr_for_sno_grouped,
                                                    start=1):
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
# pylint: enable=invalid-name

def _for_sort(x):
    return -x[1], x[0]

def plot_signature_logo(descr_full, num_res_considered=4, title=None):
    descr_filtered = descr_full.filter(['relative_sno', 'res'])
    min_sno = min(descr_filtered.relative_sno)
    descr_sno_grouped = descr_filtered.groupby('relative_sno')
    to_logo = [[] for _ in range(len(descr_sno_grouped))]
    for i, (sno, descr_per_sno) in enumerate(descr_sno_grouped):
        del sno
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
    seq_logo.Logo(to_logo, min_sno, title=title)

# pylint: disable=invalid-name, dangerous-default-value
def plot_signature_bar(descr, num_res_considered=2,
                       AA3_to_AA1=generic.AA3_to_AA1):
    """
    :param num_res_considered: Number of top-count res to show per relative
    sno position.
    """
    descr = descr.filter(['relative_sno', 'res'])
    plt.figure()
    plt.suptitle("Signatures")
    ax = plt.subplot(111)

    # Generate cmap for plots
    unique_res = np.unique(descr.res)
    num_plots = len(unique_res)
    colormap = dict()
    cmap = plt.get_cmap('gist_ncar')
    for i, res in enumerate(unique_res):
        colormap[res] = cmap(i/num_plots)

    descr = descr.groupby('relative_sno')
    bottom_val = dict()
    for sno, descr_per_sno in descr:
        res_unsorted = descr_per_sno.res.values
        res_unique = list(np.unique(res_unsorted, return_counts=True))
        res_sorted = list(zip(*sorted(zip(*res_unique), key=_for_sort)))
        res_names, res_counts = res_sorted[0][:num_res_considered], \
                                res_sorted[1][:num_res_considered]
        for name, count in zip(res_names, res_counts):
            if sno in bottom_val:
                ax.bar(sno, count, bottom=bottom_val[sno], label=name,
                       color=colormap[name])
                bottom_val[sno] += count
            else:
                ax.bar(sno, count, label=name, color=colormap[name])
                bottom_val[sno] = count

    handles, labels = ax.get_legend_handles_labels()

    # Convert to single-letter code (AA3=>AA1)
    single_letter_labels = []
    for label in labels:
        single_letter_labels.append(AA3_to_AA1[label] + " / " + label)
    labels = single_letter_labels

    # Sort labels
    labels, handles = zip(*sorted(zip(labels, handles), key=itemgetter(0)))
    labels = list(labels)
    handles = list(handles)

    # Remove duplicates
    history = []
    assert len(labels) == len(handles)
    for i in range(len(labels))[::-1]:
        if labels[i] not in history:
            history.append(labels[i])
        else:
            del labels[i]
            del handles[i]
    ax.legend(handles, labels)
# pylint: enable=invalid-name, dangerous-default-value

def plot_dihedral(descr_full, remove_labels=True, add_filename=False):
    descr_grouped = descr_full.groupby('relative_sno')

    num_plots = len(descr_grouped) - 2
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)

    plt.figure()
    plt.suptitle("Dihedrals across different relative sno position")
    for ax_count, (relative_sno, df_per_sno) in enumerate(descr_grouped):
        if ax_count in (0, len(descr_grouped) - 1):
            # Screen off first, last position, invalid dihedral values
            continue
        phis = df_per_sno.phi.values
        psis = df_per_sno.psi.values
        filenames = df_per_sno.filename.values
        seq_markers = df_per_sno.seq_marker.values

        ax = plt.subplot(num_row, num_col, ax_count)
        ax.set_title(str(relative_sno))
        ax.scatter(phis, psis, marker='x', s=1)
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])

        if remove_labels:
            ax.set_xlabel("")
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

        if add_filename:
            for phi, psi, filename, seq_marker \
                in zip(phis, psis, filenames, seq_markers):
                filename_str = str(filename) + " " + str(seq_marker)
                ax.text(phi, psi, filename_str)

def plot_hbonds_raw(descr):
    """
    Convert to heatmap eventually, maybe.
    """
    descr = descr.groupby('relative_sno')
    to_plots = dict()
    for relative_sno, df_per_sno in descr:
        donor = df_per_sno.donor.values
        counts = [len(i) for i in donor]
        to_plots[relative_sno] = (np.mean(counts), np.std(counts))
    relative_snos = list(to_plots.keys())
    means = list(zip(*to_plots.values()))[0]
    stderrs = list(zip(*to_plots.values()))[1]
    plt.figure()
    plt.suptitle("Num Hbonds")
    plt.errorbar(relative_snos, means, yerr=stderrs)

def plot_hbonds_percentage(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    to_plots = defaultdict(list)
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        counts = [len(i) for i in descr_per_set.donor.values]
        total = sum(counts)
        percent = [count/total*100 for count in counts]
        for i, sno in enumerate(relative_sno):
            to_plots[sno].append(percent[i])
    relative_snos = list(to_plots.keys())
    percents = list(to_plots.values())
    means = list([np.mean(i) for i in percents])

    plt.figure()
    plt.suptitle("Num Hbonds as Percent across Fragment")
    plt.scatter(relative_snos, means)

def plot_contacts(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        contacts = descr_per_set.contact.values
        relative_sno = descr_per_set.relative_sno.values
        assert len(contacts) == len(relative_sno)
        if not np.isnan(contacts[0]):
            for contact, sno in zip(contacts, relative_sno):
                contacts_by_pos[sno] += contact
    plt.figure()
    plt.suptitle("Contacts summed counts")
    plt.scatter(contacts_by_pos.keys(), contacts_by_pos.values())

def plot_covalents(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    covalents_by_pos = Counter()
    for _, descr_per_set in descr:
        covalents = descr_per_set.covalent.values
        relative_sno = descr_per_set.relative_sno.values
        assert len(covalents) == len(relative_sno)
        if not np.isnan(covalents[0]):
            for covalent, sno in zip(covalents, relative_sno):
                covalents_by_pos[sno] += covalent
    plt.figure()
    plt.suptitle("Covalent summed counts")
    plt.scatter(covalents_by_pos.keys(), covalents_by_pos.values())

def plot_h_contacts(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    h_contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        h_contacts = descr_per_set.h_contacts.values
        assert len(h_contacts) == len(relative_sno)
        if not np.isnan(h_contacts[0]):
            for h_contact, sno in zip(h_contacts, relative_sno):
                h_contacts_by_pos[sno] += h_contact
    plt.figure()
    plt.suptitle("H-Contact summed counts")
    plt.scatter(h_contacts_by_pos.keys(), h_contacts_by_pos.values())

def main():
    """
    Input df should have these keys:
    ['sno', 'contact', 'covalent', 'phi', 'psi', 'region', 'ss', 'ext', 'role',
     'category', 'donor', 'acc', 'res', 'CA', 'filename', 'seq_marker', 'cid']
    """
    with open(os.path.join(paths.store_dir, "descrs.pkl"),
              "rb") as pklfile:
        df = pkl.load(pklfile)
    # df.sort_index(inplace=True)
    plot_signature_logo(df)
    # plot_signature_bar(df)
    # plot_CA(df)
    plot_dihedral(df)
    plot_hbonds_raw(df)
    # plot_hbonds_percentage(df)
    plot_contacts(df)
    plot_covalents(df)
    # plot_h_contacts(df)
    # plot_dihedral_for_diff_res(df, 0)
    plt.show()


if __name__ == '__main__':
    main()
