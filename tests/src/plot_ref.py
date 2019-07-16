import logging
import pickle

import matplotlib.pyplot as plt
import numpy as np

from tests.src import paths_test
from utils import seq_logo

def plot_all():
    with open(paths_test.REF_RUN_ALL_1, 'rb') as file:
        seq_motif_map = pickle.load(file)
    plot_labelled_motifs(seq_motif_map, paths_test.REF_CREATE_SEQ_1,
                         title="Prosite_meme")
    with open(paths_test.REF_RUN_ALL_2, 'rb') as file:
        seq_motif_map = pickle.load(file)
    plot_labelled_motifs(seq_motif_map, paths_test.REF_CREATE_SEQ_1,
                         title="Prosite_mast")
    with open(paths_test.REF_RUN_ALL_3, 'rb') as file:
        seq_motif_map = pickle.load(file)
    plot_labelled_motifs(seq_motif_map, paths_test.REF_CREATE_SEQ_2,
                         title="Ioncom_mast")


def plot_labelled_motifs(seq_motif_map, fasta_filename, title=None):
    """
    This debugging function displays the sequence logo of the identified
    motifs.
    """
    motifs = _get_labelled_motifs(seq_motif_map, fasta_filename)
    converted_motifs = []
    extracted_motifs = np.array([list(i) for i in motifs])
    for AA_per_pos in extracted_motifs.T:
        AA, counts = np.unique(AA_per_pos, return_counts=True)
        sorted_i = np.argsort(counts)
        counts = counts[sorted_i]
        AA = AA[sorted_i]
        probs = counts / np.sum(counts)
        motif_this_pos = []
        for aa, prob in zip(AA, probs):
            motif_this_pos.append([aa, prob])
        converted_motifs.append(motif_this_pos)
    seq_logo.Logo(converted_motifs, -1, convert_AA3=False, title=title)


def _get_labelled_motifs(seq_motif_map, fasta_filename):
    labelled_motifs = []
    pname = None
    with open(fasta_filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                pname = line[1:5]
                continue
            if pname in seq_motif_map:
                motif_pos = seq_motif_map[pname]['sno_markers']
            else:
                #   Might have been dropped because of gapped motif, etc
                pname = None
                continue
            for pos in motif_pos:
                if len(line) <= pos + 13:
                    logging.error(
                        f"Fasta seq is shorter than pos+13, for pos in "
                        f"motif_pos. Fasta_seq: <{line}>, "
                        f"motif_pos: <{motif_pos}>, illegal pos: <{pos}>.")
                    raise Exception
                labelled_motifs.append(line[pos - 1:pos + 12])
    return labelled_motifs

if __name__ == '__main__':
    plot_all()
    plt.show()