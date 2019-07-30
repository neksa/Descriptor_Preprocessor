import find_cid_from_pname
from utils import logs

from config import paths
import preprocess
import os
from collections import Counter, OrderedDict

# todo: 1. From this: https://prosite.expasy.org/PS00164#TP, Get the aligned
#  sequences (the short pattern) in .fasta file, run meme on it to construct
#  a meme.txt, use mast on the meme txt, cut out the left-and-right portion
#  to make a length-30, and use that as my motif_pos.

# todo: meme takes forever. An alternative is to convert the aligned seqs
#  into a matrix, then call matrix2meme.

# todo: need to build a composition file using the full seq first.

import numpy as np

from utils import generic, build_composition, build_meme_from_aligned, \
    get_pname_seq

# def derive_motif_pos_from_prosite():
#     seq_file = "something"
#     aligned_seq_file = "something"
#
#     composition_file = "something"
#     build_composition.build(seq_file, composition_file)
#
#     meme_output = "meme.txt"
#     build_meme_from_aligned.build(aligned_seq_file, composition_file,
#                                  meme_output)
#
#     pname_seq_map = get_pname_seq.parse(seq_file)
#     download_pdb_files.download(pname_seq_map)
#     download_pdb_files.trim_pname_cid(pname_seq_map, paths.PDB_FOLDER)
#     pname_cid_map = find_cid_from_pname.find(pname_seq_map, paths.PDB_FOLDER)
#     with open(pname_cid_map, "wb") as file:
#
#
#     motif_pos_path = "something"
#
#     find_motifs_mast(pname_cid_path, seq_file, meme_output, 30, motif_pos_path, meme_folder=paths.MEME_MAST_FOLDER)
#     with open(motif_pos_path, 'rb') as file:
#       print(pickle.load(file))
#
#     print("\n")
#     load_pdb_info(motif_pos_path, output)
#     if storage_path is None:
#         shutil.move(conv_meme_file, paths.TRASH)
#         shutil.move(meme_folder, paths.TRASH)
#         shutil.move(motif_pos_path, paths.TRASH)


def matrix_builder(aligned_path):
    alphabets = set(generic.AA3_to_AA1.values())
    AA_to_index = {AA: i for i, AA in enumerate(sorted(alphabets))}
    matrix_counter = None
    with open(aligned_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                continue
            line = line.strip().upper()
            if matrix_counter is None:
                matrix_counter = np.zeros((len(line), len(alphabets)),
                                          dtype=int)
            for i, char in enumerate(line):
                try:
                    AA_index = AA_to_index[char]
                except KeyError:

                    # Key not found for some reason
                    continue
                matrix_counter[i, AA_index] += 1

    # matrix_counter = None
    # with open(aligned_path, 'r') as file:
    #     for line in file:
    #         if line.startswith(">"):
    #             continue
    #         line = line.strip().upper()
    #         if matrix_counter is None:
    #             matrix_counter = [Counter() for __ in range(len(line))]
    #         for i, char in enumerate(line):
    #             matrix_counter[i].update([char])
    # matrix_ordered = []
    # for counter in matrix_counter:
    #     matrix_ordered.append(OrderedDict(sorted(counter.items())))
    return matrix_counter

def write_matrix_file(matrix_ordered, output):
    """
    For meme_suite matrix2meme
    """
    output_lines = []
    for AA_counts in matrix_ordered:
        output_lines.append(" ".join(str(i) for i in AA_counts))
    single_str_line = "\n".join(output_lines)
    generic.warn_if_exist(output)
    with open(output, 'w') as file:
        file.write(single_str_line)



import pickle

def main():
    with open(os.path.join(paths.ROOT, "motif_pos.pkl"), 'rb') as file:
        print(len(pickle.load(file)))
    # preprocess.run_prosite_aligned(paths.PROSITE_ENOLASE_SEQS,
    #                                paths.PROSITE_ALIGNED_SEQS,
    #                                os.path.join(paths.ROOT,
    #                                             "output_motif_pos_new.pkl"))
    # preprocess.run_prosite_mast(paths.PROSITE_EXTRACT, 30,
    #                             paths.REF_MEME_TXT,
    #                  os.path.join(paths.ROOT, "prosite_mast_motif_pos.pkl"))
    # preprocess.run_prosite_aligned(paths.PROSITE_ENOLASE_SEQS,
    #                                paths.PROSITE_ALIGNED_SEQS,
    #                                os.path.join(paths.ROOT,
    #                                             "output_motif_pos.pkl"))

    # seq_path = os.path.join(paths.USER_INPUT, "aligner.txt")
    # matrix_ordered = matrix_builder(seq_path)
    # write_matrix_file(matrix_ordered,
    #                   os.path.join(paths.ROOT, "mine__.txt"))
    # print(matrix_ordered)
    # import sys
    # sys.exit()
    # logs.set_logging_level()
    #
    # from utils import get_pname_seq
    #
    # # pname_seq_map = get_pname_seq.parse(seq_path)
    # from meme_suite import meme_interface
    # meme_interface.create_meme_from_aligned(seq_path, 14, paths.MEME_MAST_FOLDER, num_p=7)
    # meme_txt = os.path.join(paths.MEME_MAST_FOLDER, "meme.txt")


    # logs.set_logging_level()
    # extract_path = paths.PROSITE_EXTRACT
    # motif_len = 13
    # output = paths.PID_PDB_MAP
    # num_p = 7
    # preprocess.run_prosite_meme(extract_path, motif_len, output, num_p)
    #
    # extract_path = paths.PROSITE_EXTRACT
    # motif_len = 13
    # ref_meme_txt = paths.REF_MEME_TXT
    # output = paths.PID_PDB_MAP
    # preprocess.run_prosite_mast(extract_path, motif_len, ref_meme_txt, output)
    #
    # extract_path = paths.IONCOM_EXTRACT
    # motif_len = 13
    # ref_meme_txt = paths.REF_MEME_TXT
    # output = paths.PID_PDB_MAP
    # preprocess.run_ioncom_mast(extract_path, motif_len, ref_meme_txt, output)
    pass


if __name__ == "__main__":
    main()