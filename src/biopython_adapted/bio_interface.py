import numpy as np
# from Bio.motifs import mast

from biopython_adapted import meme, mast
from utils import generic
from config import paths
from meme_suite import meme_interface

import os


def parse_meme_file(meme_path):
    generic.quit_if_missing(meme_path)
    with open(meme_path, 'r') as handle:
        motifs = meme.read(handle)
    counts = []
    motif_length = motifs[0].length
    assert motif_length >= 1
    for motif in motifs:
        motif_counts = np.array([i for i in motif.counts.values()], dtype=int)
        assert len(motif_counts) == 20
        assert len(motif_counts[0]) == motif_length
        counts.append(motif_counts)
    return counts


def parse_mast_file(meme_path):
    generic.quit_if_missing(meme_path)
    with open(meme_path, 'r') as handle:
         motifs = mast.read(handle)
    diagrams = motifs.diagrams
    return diagrams

# from tests.src import paths_test
# debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
# input_seqfile = paths_test.MEME_TEST_SEQ
# input_memefile = paths_test.REF_MEME_TXT
# output = os.path.join(debug_folder, "output_meme")
#
# meme_interface.run_mast(input_memefile, input_seqfile,
#                         output)
# mast_txt_path = os.path.join(output, "mast.txt")
# act_diagrams = parse_mast_file(mast_txt_path)
#
# from tests.src import paths_test
# meme_txt = os.path.join(paths_test.DEBUG, '2019-07-23T12:47:45.516219',
#                         'output_meme', 'meme.txt')
# fasta_filename = os.path.join(paths_test.USER_INPUT, 'meme_test_seq.fasta')
# mast_output = os.path.join(paths.ROOT, 'tmp')
# # os.mkdir(mast_output)
# meme_interface.run_mast(meme_txt, fasta_filename, mast_output)
# parse_mast_file(os.path.join(paths.INTERNAL, 'meme_mast', 'mast.txt'))
