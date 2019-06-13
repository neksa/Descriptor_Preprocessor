from collections import defaultdict, Counter
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from operator import itemgetter
import pickle as pkl
import seaborn
import os


# todo: settle magnesium now.

_aa_index = [('ALA', 'A'),
             ('CYS', 'C'),
             ('ASP', 'D'),
             ('GLU', 'E'),
             ('PHE', 'F'),
             ('GLY', 'G'),
             ('HIS', 'H'),
             ('HSE', 'H'),
             ('HSD', 'H'),
             ('ILE', 'I'),
             ('LYS', 'K'),
             ('LEU', 'L'),
             ('MET', 'M'),
             ('MSE', 'M'),
             ('ASN', 'N'),
             ('PRO', 'P'),
             ('GLN', 'Q'),
             ('ARG', 'R'),
             ('SER', 'S'),
             ('THR', 'T'),
             ('VAL', 'V'),
             ('TRP', 'W'),
             ('TYR', 'Y')]

def plot_dihedral_for_diff_res(df, sno_position):
    df = df.groupby('relative_sno')
    group = df.get_group(sno_position)

    group = group.groupby('res')
    num_plots = len(group)
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)
    ax_count = 1
    for res, group_df in group:

        phi = group_df.phi.values
        psi = group_df.psi.values

        ax = plt.subplot(num_row, num_col, ax_count)
        ax_count += 1
        ax.set_title(str(res))
        ax.scatter(phi, psi, marker='.', s=1)
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])
    print(ax_count)



def plot_CA(df):
    df = df.groupby('filename')
    for filename, df_per_file in df:
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
    return

def _for_sort(x):
    return -x[1], x[0]

def plot_signature_logo(df, num_res_considered=4):
    from utils.seq_logo import Logo
    df = df.filter(['relative_sno', 'res'])
    min_sno = min(df.relative_sno)
    df = df.groupby('relative_sno')
    to_logo = [[] for __ in range(len(df))]
    for i, (sno, df_per_sno) in enumerate(df):
        res_unsorted = df_per_sno.res.values
        res_unique = list(np.unique(res_unsorted, return_counts=True))
        total = sum(res_unique[1])
        res_sorted = list(zip(*sorted(zip(*res_unique), key=_for_sort)))
        res_names, res_counts = res_sorted[0][:num_res_considered], \
                                res_sorted[1][:num_res_considered]
        res_names = res_names[::-1]
        res_counts = res_counts[::-1]
        res_percents = list([i/total for i in res_counts])
        for name, percent in zip(res_names, res_percents):
            to_logo[i].append((name, percent))
    Logo(to_logo, min_sno)
    return


def plot_signature_bar(df, num_res_considered=2):
    """
    :param num_res_considered: Number of top-count res to show per relative
    sno position.
    """
    df = df.filter(['relative_sno', 'res'])
    plt.figure()
    plt.suptitle("Signatures")
    ax = plt.subplot(111)

    # Generate cmap for plots
    unique_res = np.unique(df.res)
    num_plots = len(unique_res)
    colormap = dict()
    cmap = plt.get_cmap('gist_ncar')
    for i, res in enumerate(unique_res):
        colormap[res] = cmap(i/num_plots)

    df = df.groupby('relative_sno')
    bottom_val = dict()
    for sno, df_per_sno in df:
        res_unsorted = df_per_sno.res.values
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

    AA3_TO_AA1 = dict(_aa_index)
    _labels = []
    for label in labels:
        _labels.append(AA3_TO_AA1[label] + " / " + label)
    labels = _labels

    # Sort labels
    labels, handles = zip(*sorted(zip(labels, handles), key=itemgetter(0)))
    labels = list(labels)
    handles = list(handles)

    # Remove duplicates
    _history = []
    assert len(labels) == len(handles)
    for i in range(len(labels))[::-1]:
        if labels[i] not in _history:
            _history.append(labels[i])
        else:
            del labels[i]
            del handles[i]

    ax.legend(handles, labels)
    return

def plot_dihedral(df, remove_labels=True, add_filename=False):
    df = df.groupby('relative_sno')

    num_plots = len(df) - 2
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)

    plt.figure()
    plt.suptitle("Dihedrals across different relative sno position")
    for ax_count, (relative_sno, df_per_sno) in enumerate(df):
        if ax_count == 0 or ax_count == len(df) - 1:
            # Screen off first, last position, invalid dihedral values
            continue
        # if ax_count not in (10, 11):
        #     continue
        phi = df_per_sno.phi.values
        psi = df_per_sno.psi.values
        filename = df_per_sno.filename.values
        seq_marker = df_per_sno.seq_marker.values

        ax = plt.subplot(num_row, num_col, ax_count)
        ax.set_title(str(relative_sno))
        ax.scatter(phi, psi, marker='x', s=1)
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])
        assert len(phi) == len(psi) == len(filename)

        if remove_labels:
            ax.set_xlabel("")
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

        if add_filename:
            for i in range(len(phi)):
                filename_str = str(filename[i]) + " " + str(seq_marker[i])
                ax.text(phi[i], psi[i], filename_str)
    return

def plot_hbonds_raw(df):
    """
    Convert to heatmap eventually, maybe.
    """
    df = df.groupby('relative_sno')
    to_plots = dict()
    for relative_sno, df_per_sno in df:
        donor = df_per_sno.donor.values
        counts = [len(i) for i in donor]
        to_plots[relative_sno] = (np.mean(counts), np.std(counts))
    relative_snos = list(to_plots.keys())
    means = list(zip(*to_plots.values()))[0]
    stderrs = list(zip(*to_plots.values()))[1]

    plt.figure()
    plt.suptitle("Num Hbonds")
    plt.errorbar(relative_snos, means, yerr=stderrs)
    return

def plot_hbonds_percentage(df):
    df = df.groupby(['filename', 'cid', 'seq_marker'])
    to_plots = defaultdict(list)
    for __, df_per_set in df:
        relative_sno = df_per_set.relative_sno.values
        counts = [len(i) for i in df_per_set.donor.values]
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
    return


# def plot_contacts(df):
#     df = df.groupby(['filename', 'cid', 'seq_marker'])
#     contacts_by_pos = Counter()
#     covalents_by_pos = Counter()
#     for __, df_per_set in df:
#         contacts = df_per_set.contact.values
#         relative_sno = df_per_set.relative_sno.values
#         assert len(contacts) == len(relative_sno)
#         if not np.isnan(contacts[0]):
#             for i in range(len(contacts)):
#                 contact = contacts[i]
#                 sno = relative_sno[i]
#                 contacts_by_pos[sno] += contact
#
#         covalents = df_per_set.covalent.values
#         assert len(covalents) == len(relative_sno)
#         if not np.isnan(covalents[0]):
#             for i in range(len(covalents)):
#                 covalent = covalents[i]
#                 sno = relative_sno[i]
#                 covalents_by_pos[sno] += covalent
#
#     plt.figure()
#     plt.suptitle("Contacts summed counts")
#     plt.scatter(contacts_by_pos.keys(), contacts_by_pos.values())
#
#     plt.figure()
#     plt.suptitle("Covalent summed counts")
#     plt.scatter(covalents_by_pos.keys(), covalents_by_pos.values())


def plot_contacts(df):
    df = df.groupby(['filename', 'cid', 'seq_marker'])
    contacts_by_pos = Counter()
    covalents_by_pos = Counter()
    # h_contacts_by_pos = Counter()
    for __, df_per_set in df:
        contacts = df_per_set.contact.values
        relative_sno = df_per_set.relative_sno.values
        assert len(contacts) == len(relative_sno)
        if not np.isnan(contacts[0]):
            for i in range(len(contacts)):
                contact = contacts[i]
                sno = relative_sno[i]
                contacts_by_pos[sno] += contact

        covalents = df_per_set.covalent.values
        assert len(covalents) == len(relative_sno)
        if not np.isnan(covalents[0]):
            for i in range(len(covalents)):
                covalent = covalents[i]
                sno = relative_sno[i]
                covalents_by_pos[sno] += covalent

        # h_contacts = df_per_set.h_contacts.values
        # assert len(h_contacts) == len(relative_sno)
        # if not np.isnan(h_contacts[0]):
        #     for i in range(len(h_contacts)):
        #         covalent = h_contacts[i]
        #         sno = relative_sno[i]
        #         h_contacts_by_pos[sno] += covalent

    plt.figure()
    plt.suptitle("Contacts summed counts")
    plt.scatter(contacts_by_pos.keys(), contacts_by_pos.values())

    plt.figure()
    plt.suptitle("Covalent summed counts")
    plt.scatter(covalents_by_pos.keys(), covalents_by_pos.values())

    # plt.figure()
    # plt.suptitle("H-Contact summed counts")
    # plt.scatter(h_contacts_by_pos.keys(), h_contacts_by_pos.values())


if __name__ == '__main__':
    import pandas as pd
    from global_config import store_dir
    """
    Input df should have these keys:
    ['sno', 'contact', 'covalent', 'phi', 'psi', 'region', 'ss', 'ext', 'role',
     'category', 'donor', 'acc', 'res', 'CA', 'filename', 'seq_marker', 'cid']
    """
    # with open(store_dir + "current.json",
    #           "r") as jsonfile:
    #     df = pd.read_json(jsonfile)

    with open(os.path.join(store_dir, "current2.pkl"),
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
    # plot_dihedral_for_diff_res(df, 0)
    plt.show()


# baserules = [
#             w.SymbolColor("GSTYC", "green", "polar"),
#             w.SymbolColor("NQ", "purple", "neutral"),
#             w.SymbolColor("KRH", "blue", "basic"),
#             w.SymbolColor("DE", "red", "acidic"),
#             w.SymbolColor("PAWFLIMV", "black", "hydrophobic")
#         ]