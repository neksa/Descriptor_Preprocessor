import sys
import find_cid_from_pname
from utils import logs

from config import paths
import preprocess
import os
import re
from collections import Counter, OrderedDict

# todo: 1. From this: https://prosite.expasy.org/PS00164#TP, Get the aligned
#  sequences (the short pattern) in .fasta file, run meme on it to construct
#  a meme.txt, use mast on the meme txt, cut out the left-and-right portion
#  to make a length-30, and use that as my motif_pos.

# todo: meme takes forever. An alternative is to convert the aligned seqs
#  into a matrix, then call matrix2meme.

# todo: need to build a composition file using the full seq first.

import numpy as np
import pickle

from utils import generic, get_pname_seq, build_composition
"""
1. Take aligned profile, build matrix
2. Download uniprot...?
3. Converge profile on uniprot...?
4. Collect sequences that match
5. Download .pdb for those that match => can we therefore run converge on 
.pdb list only?
6. Somehow cluster the sequences
7. Build for each cluster, the descr plot.

Pain points:
1. How long will running conv on uniprot take?
2. Suggestion to do a c++ rewrite of converge?


Path 1:
Prosite_extract

# todo: just rewrite conv such that its output actually make sense...? It's 
kind of lame to have to parse it separately. 

# todo: going to rewrite conv first, then when coming back here, 
fix conv_to_meme and meme_2_conv, etc. 

"""

from pdb_component import pdb_interface

from utils import clean_fasta_alphabet

from meme_suite import meme_interface
import filter_seqs
from collections import defaultdict
from converge import conv_interface
from utils import build_meme_from_aligned

# todo: check converge using the current ef-hand one, see whether original
#  maxS goes to zero as well.

def extract_meme_matrix(file):
    meme_matrix = []
    start = False
    width = None
    nsites = None
    with open(file, "r") as file:
        for line in file:
            if line.startswith("letter-probability"):
                start = True
                width = int(re.search("w= ([0-9]+)", line).group(1))
                nsites = int(re.search("nsites= ([0-9]+)", line).group(1))
                continue
            if not line.strip():
                continue
            if start:
                meme_matrix.append(list(map(float, line.strip().split(" "))))
    meme_matrix = np.array(meme_matrix)
    assert meme_matrix.shape[0] == width, meme_matrix.shape
    return (nsites, meme_matrix)



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


# def matrix_builder(aligned_path):
#     alphabets = set(generic.AA3_to_AA1.values())
#     AA_to_index = {AA: i for i, AA in enumerate(sorted(alphabets))}
#     matrix_counter = None
#     with open(aligned_path, 'r') as file:
#         for line in file:
#             if line.startswith(">"):
#                 continue
#             line = line.strip().upper()
#             if matrix_counter is None:
#                 matrix_counter = np.zeros((len(line), len(alphabets)),
#                                           dtype=int)
#             for i, char in enumerate(line):
#                 try:
#                     AA_index = AA_to_index[char]
#                 except KeyError:
#
#                     # Key not found for some reason
#                     continue
#                 matrix_counter[i, AA_index] += 1
#
#     # matrix_counter = None
#     # with open(aligned_path, 'r') as file:
#     #     for line in file:
#     #         if line.startswith(">"):
#     #             continue
#     #         line = line.strip().upper()
#     #         if matrix_counter is None:
#     #             matrix_counter = [Counter() for __ in range(len(line))]
#     #         for i, char in enumerate(line):
#     #             matrix_counter[i].update([char])
#     # matrix_ordered = []
#     # for counter in matrix_counter:
#     #     matrix_ordered.append(OrderedDict(sorted(counter.items())))
#     return matrix_counter
#
# def write_matrix_file(matrix_ordered, output):
#     """
#     For meme_suite matrix2meme
#     """
#     output_lines = []
#     for AA_counts in matrix_ordered:
#         output_lines.append(" ".join(str(i) for i in AA_counts))
#     single_str_line = "\n".join(output_lines)
#     generic.warn_if_exist(output)
#     with open(output, 'w') as file:
#         file.write(single_str_line)
#
#
def convert_nbdb_matrix_to_conv_encoder(nbdb_file, num_seqs):
    matrix = [[] for __ in range(30)]
    with open(nbdb_file, 'r') as file:
        for i, line in enumerate(file):
            probs = line.split(" ")[:20]
            for prob in probs:
                matrix[i].append(round(float(prob) * num_seqs))
    for i in matrix:
        if sum(i) > num_seqs:
            if i[0] > 0:
                i[0] -= sum(i) - num_seqs
            else:
                i[1] -= sum(i) - num_seqs
        elif sum(i) < num_seqs:
            i[0] += num_seqs - sum(i)
    for i in matrix:
        assert sum(i) == num_seqs
    return matrix
#
#
#
#
def extract_subseq_from_seqs(motif, seq, init_motif_len, final_motif_len):
    left_spacing = final_motif_len // 2 - init_motif_len // 2
    right_spacing = final_motif_len // 2 + init_motif_len // 2
    if len(seq) < motif + right_spacing:
        return None
    if motif < left_spacing:
        return None
    aligned_seq = seq[motif - left_spacing:motif + right_spacing]
    assert len(aligned_seq) == final_motif_len
    return aligned_seq
#
# def parse_conv_output(conv_output):
#     """
#     Pending conv_rewrite rewrite.
#     Essentially, we only need 1. nsites, and 2. a single matrix, and 3.
#     length of matrix. (aka width).
#     alphabet is standardised across entire project so.
#     :return:
#     """
#     # conv_output is conv_output path
#     nsite = 1
#     matrix_length = 50
#     matrix = np.zeros((matrix_length, 20), dtype=float)
#     return matrix_length, nsite, matrix

#
#
# def post_meme():
#     PDB_SEQ_FILE = os.path.join(paths.ROOT, 'whatever.fasta')
#     MEME_TXT = 'whatever.txt'
#     # For uniprot, not pdb_seq. Because converge/meme generation is based on
#     # uniprot.
#     COMPOSITION_FILE = os.path.join(paths.ROOT, 'whatever.fasta')
#     MOTIF_POS_OUTPUT = os.path.join(paths.USER_OUTPUT, "motif_pos.pkl")
#     # composition_file should be provided, can be generated as below:
#
#     # Intermediate paths
#     MEME_FOLDER = paths.MEME_MAST_FOLDER
#     # can someone do the cleaning in advance?
#     CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(2)
#     PDB_SEQ_FROM_PDB_FILES = paths.TMP_FILE_TEMPLATE.format(3)
#
#     meme_interface.run_mast(MEME_TXT, PDB_SEQ_FILE, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     # this should give me both cid, and pdb_id. It doesn't now.
#     pdb_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
#                                                        matrix_length)
#     pdb_seqs = dict()
#     for pdb_id, (cid, _motif) in pdb_motif_map.items():
#         seq = pdb_interface.get_seq_for(pdb_id, cid)
#         pdb_seqs[(pdb_id, cid)] = seq
#     write_as_seq(pdb_seqs, PDB_SEQ_FROM_PDB_FILES)
#     meme_interface.run_mast(MEME_TXT, PDB_SEQ_FROM_PDB_FILES, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     # this should give me both cid, and pdb_id. It doesn't now.
#     pdb_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
#                                                        matrix_length)
#     output_motif_pos(pdb_motif_map, MOTIF_POS_OUTPUT)
#
# def post_conv():
#     """
#     I need to take converge output,
#     1. convert it to meme format
#     2. run it on pdb_seqs
#     3. Collect all pdb_id, cid that match
#     4. Download relevant pdb files
#     5. Rebuild another seq_file based on collected pdb files, and cid.
#     6. Run mast again
#     7. Collect the correct motif_pos, cid, pdb_id
#     8. Output as motif_pos.pkl, for descr_calculator to take over
#
#     Reason for this is:
#     Order of sequences in pdb_seqs.fasta may not correspond to that in .pdb
#     files. So, running mast on the recompiled seq_file is safest.
#     """
#     PDB_SEQ_FILE = os.path.join(paths.ROOT, 'whatever.fasta')
#     # For uniprot, not pdb_seq. Because converge/meme generation is based on
#     # uniprot.
#     COMPOSITION_FILE = os.path.join(paths.ROOT, 'whatever.fasta')
#     MOTIF_POS_OUTPUT = os.path.join(paths.USER_OUTPUT, "motif_pos.pkl")
#     # composition_file should be provided, can be generated as below:
#
#     # Intermediate paths
#     MEME_FOLDER = paths.MEME_MAST_FOLDER
#     MEME_TXT = paths.TMP_FILE_TEMPLATE.format(1)
#     # can someone do the cleaning in advance?
#     CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(2)
#     PDB_SEQ_FROM_PDB_FILES = paths.TMP_FILE_TEMPLATE.format(3)
#
#     matrix_length, nsite, matrix = parse_conv_output(CONV_OUTPUT)
#     write_to_meme(MEME_TXT, matrix_length, nsite, matrix)
#     meme_interface.run_mast(MEME_TXT, PDB_SEQ_FILE, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     # this should give me both cid, and pdb_id. It doesn't now.
#     pdb_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
#                                                        matrix_length)
#     pdb_seqs = dict()
#     for pdb_id, (cid, _motif) in pdb_motif_map.items():
#         seq = pdb_interface.get_seq_for(pdb_id, cid)
#         pdb_seqs[(pdb_id, cid)] = seq
#     write_as_seq(pdb_seqs, PDB_SEQ_FROM_PDB_FILES)
#     meme_interface.run_mast(MEME_TXT, PDB_SEQ_FROM_PDB_FILES, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     # this should give me both cid, and pdb_id. It doesn't now.
#     pdb_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
#                                                        matrix_length)
#     output_motif_pos(pdb_motif_map, MOTIF_POS_OUTPUT)
#
#
#
# def seq_aligned_meme():
#     # aligned sequences need to have no gap
#     # todo: test for that?
#
#     # Known matched sequences
#     SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_seqs.fasta')
#     ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_aligned.fasta')
#     # Length of actual motif
#     FINAL_MOTIF_LEN = 50
#     OUTPUT_MEME = os.path.join(paths.USER_OUTPUT, "seq_aligned_meme.txt")
#     # composition_file should be provided, can be generated as below:
#     COMPOSITION_FILE = os.path.join(paths.ROOT, "uniprot_full_composition.txt")
#     build_composition.build("uniprot_full.fasta", COMPOSITION_FILE)
#
#     # Intermediate paths
#     MEME_FOLDER = paths.MEME_MAST_FOLDER
#     MEME_TXT = paths.TMP_FILE_TEMPLATE.format(1)
#     CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(2)
#     TMP_ALIGNED_SEQS = paths.TMP_FILE_TEMPLATE.format(3)
#
#     INIT_MOTIF_LEN = preprocess._get_motif_len_from_aligned(ALIGNED_SEQ_FILE)
#     clean_fasta_alphabet.screen(fasta_path=SEQ_FILE, output_path=CLEANED_SEQ)
#     filter_seqs.delete_short_seqs(CLEANED_SEQ, threshold=FINAL_MOTIF_LEN)
#
#     build_meme_from_aligned.build(ALIGNED_SEQ_FILE, MEME_TXT, COMPOSITION_FILE)
#
#     meme_interface.run_mast(MEME_TXT, CLEANED_SEQ, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
#                                                                INIT_MOTIF_LEN)
#
#     acc_seq_map = get_pname_seq.parse_raw(CLEANED_SEQ)
#     aligned_sequences = []
#     for acc, motifs in acc_motif_map.items():
#         if acc not in acc_seq_map:
#             continue
#         for motif in motifs:
#             subseq = extract_subseq_from_seqs(motif, acc_seq_map[acc],
#                                               INIT_MOTIF_LEN, FINAL_MOTIF_LEN)
#             if subseq:
#                 aligned_sequences.append(subseq)
#
#     with open(TMP_ALIGNED_SEQS, "w") as file:
#         for seq in aligned_sequences:
#             file.write(">RAND\n")
#             file.write(seq + "\n")
#
#     build_meme_from_aligned.build(TMP_ALIGNED_SEQS, OUTPUT_MEME, COMPOSITION_FILE)
#
#     # end
#     # to check logo of this, run ./meme/ceqlogo -i1 output_meme_mine.txt -o
#     # logo.eps -f EPS
#     # output as matrix, and number of Kmatches.
#     # Number of Kmatches is approximately number of aligned seqs
#
#
#


def seq_aligned_converge_testing():
    '''
    16/10/2019
    Testing to see if my conv_rewrite gives same output as the original
    converge. Mainly in tracking of maxS.

    To do this, we'll take aligned_sequences, run make_meme on it, then run
    meme_to_conv, then use that conv as seed_matrix, and run conv_orig on the
    full uniprot data.
    :return:
    '''
    # aligned sequences need to have no gap
    # todo: test for that?

    # Known matched sequences
    SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_seqs.fasta')
    ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_aligned.txt')
    # Length of actual motif
    FINAL_MOTIF_LEN = 30
    OUTPUT_MATRIX = os.path.join(paths.USER_OUTPUT,
                                 "seq_aligned_converge_output.txt")
    # composition_file should be provided, can be generated as below:
    COMPOSITION_FILE = os.path.join(paths.ROOT, "uniprot_full_composition.txt")
    build_composition.build(os.path.join(paths.ROOT, "uniprot_full.fasta"), \
                            COMPOSITION_FILE)

    # Intermediate paths
    MEME_FOLDER = paths.MEME_MAST_FOLDER
    MEME_TXT = paths.TMP_FILE_TEMPLATE.format(1)
    CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(2)

    INIT_MOTIF_LEN = preprocess._get_motif_len_from_aligned(ALIGNED_SEQ_FILE)
    clean_fasta_alphabet.screen(fasta_path=SEQ_FILE, output_path=CLEANED_SEQ)
    filter_seqs.delete_short_seqs(CLEANED_SEQ, threshold=FINAL_MOTIF_LEN)

    build_meme_from_aligned.build(ALIGNED_SEQ_FILE, MEME_TXT, COMPOSITION_FILE)

    meme_interface.run_mast(MEME_TXT, CLEANED_SEQ, MEME_FOLDER)
    mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
    acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
                                                               INIT_MOTIF_LEN)

    acc_seq_map = get_pname_seq.parse_raw(CLEANED_SEQ)
    aligned_sequences = []
    for acc, motifs in acc_motif_map.items():
        if acc not in acc_seq_map:
            continue
        for motif in motifs:
            subseq = extract_subseq_from_seqs(motif, acc_seq_map[acc],
                                              INIT_MOTIF_LEN, FINAL_MOTIF_LEN)
            if subseq:
                aligned_sequences.append(subseq)


    composition_file = os.path.join(paths.ROOT, 'uniprot_full_composition.txt')
    aligned_seq_file = os.path.join(paths.ROOT, "aligned")
    with open(aligned_seq_file, 'w') as file:
        for i, seq in enumerate(aligned_sequences):
            file.write(">RAND_{}\n".format(i))
            file.write(seq + "\n")

    output_meme_for_testing = os.path.join(paths.ROOT,
                                           "output_meme_for_testing")
    build_meme_from_aligned.build(aligned_seq_file, output_meme_for_testing,
                                  composition_file)
    # sys.exit()
    # from converge.conv_interface import convert_meme_to_conv
    from converge import conv_interface
    conv_seed_mat = os.path.join(paths.ROOT, "conv_seed")
    conv_interface.convert_meme_to_conv(output_meme_for_testing,
                                        composition_file, conv_seed_mat,
                                        minimal=True)

    output_matrix = np.zeros((FINAL_MOTIF_LEN, 20), dtype=int)
    for seq in aligned_sequences:
        for i, char in enumerate(seq):
            if char in generic.AA1_TO_I:
                output_matrix[i][generic.AA1_TO_I[char]] += 1

    with open(OUTPUT_MATRIX, 'w') as file:
        for line in output_matrix:
            for value in line:
                file.write(str(value) + " ")
            file.write("\n")
    print(f"Kmatches is {len(aligned_sequences)}.")

    #     next, we run converge.  #     assumption that proteome, blosum is
    #     already encoded
    from converge import conv_interface

    OUTPUT_MEME = os.path.join(paths.ROOT, "output_meme.txt")

    # conv_interface.run(OUTPUT_MATRIX, FINAL_MOTIF_LEN, len(aligned_sequences),
    #                    OUTPUT_MEME)

    conv_interface.run(OUTPUT_MATRIX, FINAL_MOTIF_LEN, len(aligned_sequences),
                       OUTPUT_MEME)

# from converge import conv_interface
#
# OUTPUT_MEME = os.path.join(paths.ROOT, "output_meme.txt")
# composition_unused = paths.TMP_FILE_TEMPLATE.format(1)
# conv_output = paths.TMP_FILE_TEMPLATE.format("conv_output")
#
# conv_interface.convert_meme_to_conv(OUTPUT_MEME, composition_unused, conv_output,
#                                     minimal=True)
# import sys
# sys.exit()

import traceback
import logging


def seq_aligned_converge(num_seqs=None):
    """
    # todo 26092019: if this works, repeat for GxGGxG and GxxGxG, then build all
    descrs together.

    Will need to run this for GxGxxG.
    We are supplied a matrix. Will use it to produce seed_matrix to run on
    converge on uniprot, then the rest are the same. Converged profile then
    run mast on rcsb seq file, etc.
    :return:
    """

    # Known matched sequences
    # SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_seqs.fasta')
    # ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_aligned.txt')

    # SEQ_FILE = os.path.join(paths.ROOT, 'efhand_aligned_seqs.fasta')
    # ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'efhand_aligned.txt')
    # Length of actual motif
    FINAL_MOTIF_LEN = 30
    OUTPUT_MATRIX = os.path.join(paths.USER_OUTPUT,
                                 "seq_aligned_converge_output.txt")
    # composition_file should be provided, can be generated as below:
    COMPOSITION_FILE = os.path.join(paths.ROOT, "uniprot_full_composition.txt")
    build_composition.build(os.path.join(paths.ROOT, "uniprot_full.fasta"), \
                            COMPOSITION_FILE)

    # Intermediate paths
    MEME_FOLDER = paths.MEME_MAST_FOLDER
    MEME_TXT = paths.TMP_FILE_TEMPLATE.format(1)
    CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(2)

    # INIT_MOTIF_LEN = preprocess._get_motif_len_from_aligned(ALIGNED_SEQ_FILE)
    # clean_fasta_alphabet.screen(fasta_path=SEQ_FILE, output_path=CLEANED_SEQ)
    # filter_seqs.delete_short_seqs(CLEANED_SEQ, threshold=FINAL_MOTIF_LEN)
    #
    # build_meme_from_aligned.build(ALIGNED_SEQ_FILE, MEME_TXT, COMPOSITION_FILE)

    # meme_interface.run_mast(MEME_TXT, CLEANED_SEQ, MEME_FOLDER)
    # mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
    # acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
    #                                                            INIT_MOTIF_LEN)
    #
    # acc_seq_map = get_pname_seq.parse_raw(CLEANED_SEQ)
    # aligned_sequences = []
    # for acc, motifs in acc_motif_map.items():
    #     if acc not in acc_seq_map:
    #         continue
    #     for motif in motifs:
    #         subseq = extract_subseq_from_seqs(motif, acc_seq_map[acc],
    #                                           INIT_MOTIF_LEN, FINAL_MOTIF_LEN)
    #         if subseq:
    #             aligned_sequences.append(subseq)
    #
    # output_matrix = np.zeros((FINAL_MOTIF_LEN, 20), dtype=int)
    # for seq in aligned_sequences:
    #     for i, char in enumerate(seq):
    #         if char in generic.AA1_TO_I:
    #             output_matrix[i][generic.AA1_TO_I[char]] += 1

    # This is for running using the merged pssm
    merged_pssm = os.path.join(paths.ROOT, 'written_meme.txt')
    num_seqs, output_matrix = extract_meme_matrix(merged_pssm)

    # nbdb_file = os.path.join(paths.ROOT, "GxGxxG_pssm.txt")
    # num_seqs = 4000
    # output_matrix = convert_nbdb_matrix_to_conv_encoder(nbdb_file, num_seqs)

    with open(OUTPUT_MATRIX, 'w') as file:
        for line in output_matrix:
            for value in line:
                file.write(str(value) + " ")
            file.write("\n")
    print(f"Kmatches is {num_seqs}.")

    from converge import conv_interface

    OUTPUT_MEME = os.path.join(paths.ROOT, "output_meme.txt")

    # conv_interface.run(OUTPUT_MATRIX, FINAL_MOTIF_LEN, len(aligned_sequences),
    #                    OUTPUT_MEME)
    # conv_interface.run(OUTPUT_MATRIX, FINAL_MOTIF_LEN, num_seqs,
    #                    OUTPUT_MEME)
    OUTPUT_MEME = merged_pssm

    pdb_seq_file = paths.RCSB_SEQS_FASTA
    clean_fasta_alphabet.screen(pdb_seq_file)
    meme_interface.run_mast(OUTPUT_MEME, pdb_seq_file, paths.MEME_MAST_FOLDER)
    mast_txt_path = os.path.join(paths.MEME_MAST_FOLDER, 'mast.txt')

    pdb_cid_motif_raw = _get_motif_diagram_mast_pdb(mast_txt_path)

    pdb_cid_seq = dict()
    # three_pdb = ["1b8c", "1ig5", "3tuy"]
    # three_cid = ['A', "A", "B"]

    for i, (pdb_id, cid) in enumerate(pdb_cid_motif_raw):
    # for i, (pdb_id, cid) in enumerate(zip(three_pdb, three_cid)):
        print(i)
        try:
            seq = pdb_interface.get_seq_for(pdb_id, cid=cid)
            if seq is None:
                continue
        except Exception as e:
            logging.error(f"get_seq_for() fails for pdb_id/cid {pdb_id}/{cid}. "
                          f"Skipping.")
            logging.error(f"Traceback: <{traceback.format_exc()}>")
            logging.error(f"Error_msg: <{e}>\n\n")
            continue
        pdb_cid_seq[(pdb_id, cid)] = seq
    print(pdb_cid_seq)

    with open(paths.TMP_FILE_TEMPLATE.format("tmp_pdb_cid_seq.pkl"), 'wb') as \
            file:
        pickle.dump(pdb_cid_seq, file, -1)
    pdb_cid_seq = OrderedDict(sorted(pdb_cid_seq.items()))
    pdb_seq_direct_from_pdb_files = paths.TMP_FILE_TEMPLATE.format(12)
    with open(pdb_seq_direct_from_pdb_files, 'w') as file:
        for (pdb_id, cid), seq in pdb_cid_seq.items():
            file.write(f">{pdb_id}_{cid}\n")
            file.write(seq + "\n")

    clean_fasta_alphabet.screen(pdb_seq_direct_from_pdb_files)
    filter_seqs.delete_short_seqs(pdb_seq_direct_from_pdb_files,
                                  FINAL_MOTIF_LEN)
    meme_interface.run_mast(OUTPUT_MEME, pdb_seq_direct_from_pdb_files,
                            paths.MEME_MAST_FOLDER)
    mast_txt_path = os.path.join(paths.MEME_MAST_FOLDER, 'mast.txt')

    pdb_cid_motif_raw = _get_motif_diagram_mast_pdb(mast_txt_path)
    print("pdb_cid_motif_raw")
    print(pdb_cid_motif_raw)
    pdb_cid_motif_map = meme_interface._adjust_motif_diagram(pdb_cid_motif_raw,
                                                             FINAL_MOTIF_LEN)
    # print("pdb_cid_motif_map")
    # print(pdb_cid_motif_map)
    # pdb_cid_motif_map = motif_finder._delete_gapped_motifs(pdb_cid_motif_map,
    #                                                    pdb_seq_direct_from_pdb_files)
    print("pdb_cid_motif_map2")
    print(pdb_cid_motif_map)

    motif_positions = defaultdict(dict)
    for (pdb_id, cid), motif_pos in pdb_cid_motif_map.items():
        motif_positions[pdb_id]['sno_markers'] = motif_pos
        motif_positions[pdb_id]['cid'] = cid
    motif_positions = OrderedDict(sorted(motif_positions.items()))
    print(motif_positions)
    with open(os.path.join(paths.ROOT, "motif_pos_tmp.pkl"), 'wb') as file:
        pickle.dump(motif_positions, file, -1)
    return motif_positions


#
# def seq_aligned_converge():
"""
26092019: Frozen, for running of MG and EF-hand. 
"""
#     # aligned sequences need to have no gap
#     # todo: test for that?
#
#     # Known matched sequences
#     SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_seqs.fasta')
#     ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'MG_dxdxxd_aligned.txt')
#
#     # SEQ_FILE = os.path.join(paths.ROOT, 'efhand_aligned_seqs.fasta')
#     # ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'efhand_aligned.txt')
#     # Length of actual motif
#     FINAL_MOTIF_LEN = 30
#     OUTPUT_MATRIX = os.path.join(paths.USER_OUTPUT,
#                                  "seq_aligned_converge_output.txt")
#     # composition_file should be provided, can be generated as below:
#     COMPOSITION_FILE = os.path.join(paths.ROOT, "uniprot_full_composition.txt")
#     build_composition.build(os.path.join(paths.ROOT, "uniprot_full.fasta"), \
#                                           COMPOSITION_FILE)
#
#     # Intermediate paths
#     MEME_FOLDER = paths.MEME_MAST_FOLDER
#     MEME_TXT = paths.TMP_FILE_TEMPLATE.format(1)
#     CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(2)
#
#     INIT_MOTIF_LEN = preprocess._get_motif_len_from_aligned(ALIGNED_SEQ_FILE)
#     clean_fasta_alphabet.screen(fasta_path=SEQ_FILE, output_path=CLEANED_SEQ)
#     filter_seqs.delete_short_seqs(CLEANED_SEQ, threshold=FINAL_MOTIF_LEN)
#
#     build_meme_from_aligned.build(ALIGNED_SEQ_FILE, MEME_TXT, COMPOSITION_FILE)
#
#     meme_interface.run_mast(MEME_TXT, CLEANED_SEQ, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
#                                                                INIT_MOTIF_LEN)
#
#     acc_seq_map = get_pname_seq.parse_raw(CLEANED_SEQ)
#     aligned_sequences = []
#     for acc, motifs in acc_motif_map.items():
#         if acc not in acc_seq_map:
#             continue
#         for motif in motifs:
#             subseq = extract_subseq_from_seqs(motif, acc_seq_map[acc],
#                                               INIT_MOTIF_LEN, FINAL_MOTIF_LEN)
#             if subseq:
#                 aligned_sequences.append(subseq)
#
#     output_matrix = np.zeros((FINAL_MOTIF_LEN, 20), dtype=int)
#     for seq in aligned_sequences:
#         for i, char in enumerate(seq):
#             if char in generic.AA1_TO_I:
#                 output_matrix[i][generic.AA1_TO_I[char]] += 1
#
#     with open(OUTPUT_MATRIX, 'w') as file:
#         for line in output_matrix:
#             for value in line:
#                 file.write(str(value) + " ")
#             file.write("\n")
#     print(f"Kmatches is {len(aligned_sequences)}.")
#
#     #     next, we run converge.  #     assumption that proteome, blosum is
#     #     already encoded
#     from converge import conv_interface
#
#     OUTPUT_MEME = os.path.join(paths.ROOT, "output_meme.txt")
#
#     # todo: need to take aligned_seqs, convert it to meme.txt, then to
#     #  converge matrix format, then compare between conv_orig and
#     #  converge_calculator, tracking maxS to see differences.
#
#     conv_interface.run(OUTPUT_MATRIX, FINAL_MOTIF_LEN, len(aligned_sequences),
#                        OUTPUT_MEME)
#
#     # conv_interface.run(OUTPUT_MATRIX, FINAL_MOTIF_LEN, 3522,
#     #                    OUTPUT_MEME)
#
#     # from converge import meme_to_conv
#     # matrices = []
#     # with open(OUTPUT_MATRIX, 'r') as file:
#     #     for line in file:
#     #         matrix_line = []
#     #         values = line.strip().split(" ")
#     #         for value in values:
#     #             matrix_line.append(int(value)/3522)
#     #         matrices.append(matrix_line)
#     # print(len(matrices))
#     # TMP_CONV_INPUT_MAT = paths.TMP_FILE_TEMPLATE.format("tmp_conv_input_mat")
#     #
#     # meme_to_conv._write_matrix(matrices, 3522, TMP_CONV_INPUT_MAT)
#
#     # after converge gives converged profile:
#     # 1. run mast on pdb_seq, extract list of pdb_id and cid
#     # 2. Extract seq of these pdb_id and cid from pdb files
#     # 2.5: Check for short seqs?
#     # 3. Run mast again on this set, to get motif_pos
#     # 4. Check for gapped motifs?
#     # 4. Then, output motif_pos.
#
#     pdb_seq_file = paths.RCSB_SEQS_FASTA
#     clean_fasta_alphabet.screen(pdb_seq_file)
#     meme_interface.run_mast(OUTPUT_MEME, pdb_seq_file, paths.MEME_MAST_FOLDER)
#     mast_txt_path = os.path.join(paths.MEME_MAST_FOLDER, 'mast.txt')
#
#     pdb_cid_motif_raw = _get_motif_diagram_mast_pdb(mast_txt_path)
#
#     pdb_cid_seq = dict()
#     print(len(pdb_cid_motif_raw))
#     for i, (pdb_id, cid) in enumerate(pdb_cid_motif_raw.keys()):
#         print(i)
#         try:
#             seq = pdb_interface.get_seq_for(pdb_id, cid=cid)
#             if seq is None:
#                 continue
#         except Exception as e:
#             logging.error(f"get_seq_for() fails for pdb_id/cid {pdb_id}/{cid}. "
#                           f"Skipping.")
#             logging.error(f"Traceback: <{traceback.format_exc()}>")
#             logging.error(f"Error_msg: <{e}>\n\n")
#             continue
#         pdb_cid_seq[(pdb_id, cid)] = seq
#     with open(paths.TMP_FILE_TEMPLATE.format("tmp_pdb_cid_seq.pkl"), 'wb') as\
#             file:
#         pickle.dump(pdb_cid_seq, file, -1)
#     pdb_cid_seq = OrderedDict(sorted(pdb_cid_seq.items()))
#     pdb_seq_direct_from_pdb_files = paths.TMP_FILE_TEMPLATE.format(12)
#     with open(pdb_seq_direct_from_pdb_files, 'w') as file:
#         for (pdb_id, cid), seq in pdb_cid_seq.items():
#             file.write(f">{pdb_id}_{cid}\n")
#             file.write(seq + "\n")
#     clean_fasta_alphabet.screen(pdb_seq_file)
#     filter_seqs.delete_short_seqs(pdb_seq_direct_from_pdb_files, FINAL_MOTIF_LEN)
#     meme_interface.run_mast(OUTPUT_MEME, pdb_seq_direct_from_pdb_files, paths.MEME_MAST_FOLDER)
#     mast_txt_path = os.path.join(paths.MEME_MAST_FOLDER, 'mast.txt')
#
#     pdb_cid_motif_raw = _get_motif_diagram_mast_pdb(mast_txt_path)
#     print("pdb_cid_motif_raw")
#     print(pdb_cid_motif_raw)
#     pdb_cid_motif_map = meme_interface._adjust_motif_diagram(pdb_cid_motif_raw,
#                                                                FINAL_MOTIF_LEN)
#     # print("pdb_cid_motif_map")
#     # print(pdb_cid_motif_map)
#     # pdb_cid_motif_map = motif_finder._delete_gapped_motifs(pdb_cid_motif_map,
#     #                                                    pdb_seq_direct_from_pdb_files)
#     print("pdb_cid_motif_map2")
#     print(pdb_cid_motif_map)
#
#     motif_positions = defaultdict(dict)
#     for (pdb_id, cid), motif_pos in pdb_cid_motif_map.items():
#         motif_positions[pdb_id]['sno_markers'] = motif_pos
#         motif_positions[pdb_id]['cid'] = cid
#     motif_positions = OrderedDict(sorted(motif_positions.items()))
#     print(motif_positions)
#     with open(os.path.join(paths.ROOT, "motif_pos_tmp.pkl"), 'wb') as file:
#         pickle.dump(motif_positions, file, -1)
#     return motif_positions


def _get_motif_diagram_mast_pdb(input_txt):
    """
    We obtain from meme/mast output (meme.txt/mast.txt) the motif diagram,
    showing the location of the matched motifs for each pname sequence. We
    then adjust it such that the locations are absolute, relative to the
    start of the sequence rather than from the end of the previous motif.
    Motifs with gaps in between are also deleted, mainly because the
    subsequent descriptor code assumes a continuous sequence for the len-30
    analysis segment.
    """
    pname_motif_map = dict()
    motif_diagram_sep = "[1]"
    with open(input_txt, 'r') as rfile:
        correct_segment = False
        in_area = False
        for line in rfile:
            if line.startswith("SECTION II: MOTIF DIAGRAMS"):
                correct_segment = True
                continue
            if correct_segment and line.startswith("-------------   "):
                in_area = True
                continue
            if in_area and not line.strip():
                break
            if in_area:
                pdb_id, cid, motif_diagram = line[:4], line[5], line[43:]
                if motif_diagram_sep not in motif_diagram:
                    continue
                raw_motif_positions = motif_diagram.strip().\
                    split(motif_diagram_sep)
                motif_positions = []
                for pos in raw_motif_positions[:-1]:
                    pos = pos.strip()
                    spacers = pos.replace("-", " ").strip().split(" ")
                    if len(spacers) == 0 or spacers[0] == "":
                        motif_positions.append(0)
                    else:
                        try:
                            assert len(spacers) <= 1, spacers
                            motif_positions.append(int(spacers[0]))
                        except Exception as e:
                            logging.error(
                                f"_get_motif_diagram_mast_pdb() fails for line "
                                f"{line}. Skipping.")
                            logging.error(
                                f"Traceback: <{traceback.format_exc()}>")
                            logging.error(f"Error_msg: <{e}>\n\n")
                            continue
                pname_motif_map[(pdb_id, cid)] = motif_positions
    return pname_motif_map

# def _get_motif_diagram_mast_pdb(input_txt):
#     """
#     We obtain from meme/mast output (meme.txt/mast.txt) the motif diagram,
#     showing the location of the matched motifs for each pname sequence. We
#     then adjust it such that the locations are absolute, relative to the
#     start of the sequence rather than from the end of the previous motif.
#     Motifs with gaps in between are also deleted, mainly because the
#     subsequent descriptor code assumes a continuous sequence for the len-30
#     analysis segment.
#     """
#     pname_motif_map = dict()
#     motif_diagram_sep = "-[1]-"
#     with open(input_txt, 'r') as rfile:
#         correct_segment = False
#         in_area = False
#         for line in rfile:
#             if line.startswith("SECTION II: MOTIF DIAGRAMS"):
#                 correct_segment = True
#                 continue
#             if correct_segment and line.startswith("-------------   "):
#                 in_area = True
#                 continue
#             if in_area and not line.strip():
#                 break
#             if in_area:
#                 pdb_id, cid, motif_diagram = line[:4], line[5], line[43:]
#                 if motif_diagram_sep not in motif_diagram:
#                     continue
#                 # Remove [1]_5 if first motif is right at the front.
#                 if motif_diagram.startswith("["):
#                     motif_diagram = motif_diagram[4:]
#                 raw_motif_positions = motif_diagram.split(motif_diagram_sep)
#                 motif_positions = []
#                 for pos in raw_motif_positions[:-1]:
#                     if pos == "":
#                         motif_positions.append(0)
#                     else:
#                         motif_positions.append(int(pos))
#                 pname_motif_map[(pdb_id, cid)] = motif_positions
#     return pname_motif_map
#

#
# def seq_only_converge():
#     """
#     04/09/2019
#
#     Starting from only sequences that match. Uses the first match from meme
#     as seed_motif. Align back to known sequences, expand flanks, use as
#     seed_matrix, output as conv_matrix for converge_encoder to convert to
#     serialised format for converge.
#
#     Single match from meme may not correspond to desired motif, need to
#     manually verify.
#     """
#     # Known matched sequences
#     SEQ_FILE = os.path.join(paths.ROOT, 'efhand_seq.fasta')
#     # Length of initial motif in SEQ_FILE, keep it short so match to motif is
#     # correct
#     INIT_MOTIF_LEN = 15
#     # Number of processors for initial meme run
#     NUM_P = 7
#     # Length of actual motif
#     FINAL_MOTIF_LEN = 50
#     OUTPUT_MATRIX = os.path.join(paths.USER_OUTPUT,
#                                  "seq_only_converge_output.txt")
#
#     # Intermediate paths
#     MEME_FOLDER = paths.MEME_MAST_FOLDER
#     CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(1)
#
#     # Remove extraneous symbols that would cause meme or converge_encoder
#     # to fail.
#     clean_fasta_alphabet.screen(fasta_path=SEQ_FILE, output_path=CLEANED_SEQ)
#     filter_seqs.delete_short_seqs(CLEANED_SEQ, threshold=FINAL_MOTIF_LEN)
#     meme_interface.run_meme(CLEANED_SEQ, INIT_MOTIF_LEN, MEME_FOLDER, NUM_P)
#     meme_txt = os.path.join(MEME_FOLDER, "meme.txt")
#     # This meme_txt is for finding motifs in matched_sequences
#     # To verify, run make_logo_meme().
#     # todo: build single function for making ceqlogo, add here.
#
#     # Matches meme to matched_sequences, for deriving of final motif later
#     meme_interface.run_mast(meme_txt, CLEANED_SEQ, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
#                                                                INIT_MOTIF_LEN)
#
#     acc_seq_map = get_pname_seq.parse_raw(CLEANED_SEQ)
#     aligned_sequences = []
#     for acc, motifs in acc_motif_map.items():
#         if acc not in acc_seq_map:
#             continue
#         for motif in motifs:
#             subseq = extract_subseq_from_seqs(motif, acc_seq_map[acc],
#                                               INIT_MOTIF_LEN, FINAL_MOTIF_LEN)
#             if subseq:
#                 aligned_sequences.append(subseq)
#
#     output_matrix = np.zeros((FINAL_MOTIF_LEN, 20), dtype=int)
#     for seq in aligned_sequences:
#         for i, char in enumerate(seq):
#             if char in generic.AA1_TO_I:
#                 output_matrix[i][generic.AA1_TO_I[char]] += 1
#
#     with open(OUTPUT_MATRIX, 'w') as file:
#         for line in output_matrix:
#             for value in line:
#                 file.write(str(value) + " ")
#             file.write("\n")
#     print(f"Kmatches is {len(aligned_sequences)}.")
#     return
#
#
# def seq_only_meme():
#     """
#     04/09/2019
#
#     Starting from only sequences that match. Uses the first match from meme
#     as seed_motif. Align back to known sequences, expand flanks, use as
#     seed_matrix, output as conv_matrix for converge_encoder to convert to
#     serialised format for converge.
#
#     Single match from meme may not correspond to desired motif, need to
#     manually verify.
#     """
#     # Known matched sequences
#     SEQ_FILE = os.path.join(paths.ROOT, 'efhand_seq.fasta')
#     # Length of initial motif in SEQ_FILE, keep it short so match to motif is
#     # correct
#     INIT_MOTIF_LEN = 15
#     # Number of processors for initial meme run
#     NUM_P = 7
#     # Length of actual motif
#     FINAL_MOTIF_LEN = 50
#     OUTPUT_MEME = os.path.join(paths.USER_OUTPUT,
#                                  "seq_only_meme.txt")
#     # composition_file should be provided, can be generated as below:
#     COMPOSITION_FILE = os.path.join(paths.ROOT, "uniprot_full_composition.txt")
#     build_composition.build("uniprot_full.fasta", COMPOSITION_FILE)
#
#     # Intermediate paths
#     MEME_FOLDER = paths.MEME_MAST_FOLDER
#     CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(1)
#     ALIGNED_SEQS = paths.TMP_FILE_TEMPLATE.format(2)
#
#     # Remove extraneous symbols that would cause meme or converge_encoder
#     # to fail.
#     clean_fasta_alphabet.screen(fasta_path=SEQ_FILE, output_path=CLEANED_SEQ)
#     filter_seqs.delete_short_seqs(CLEANED_SEQ, threshold=INIT_MOTIF_LEN)
#     meme_interface.run_meme(CLEANED_SEQ, INIT_MOTIF_LEN, MEME_FOLDER, NUM_P)
#     meme_txt = os.path.join(MEME_FOLDER, "meme.txt")
#     # This meme_txt is for finding motifs in matched_sequences
#     # To verify, run make_logo_meme().
#     # todo: build single function for making ceqlogo, add here.
#
#     # Matches meme to matched_sequences, for deriving of final motif later
#     meme_interface.run_mast(meme_txt, CLEANED_SEQ, MEME_FOLDER)
#     mast_txt_path = os.path.join(MEME_FOLDER, 'mast.txt')
#     acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
#                                                                INIT_MOTIF_LEN)
#
#     acc_seq_map = get_pname_seq.parse_raw(CLEANED_SEQ)
#     aligned_sequences = []
#     for acc, motifs in acc_motif_map.items():
#         if acc not in acc_seq_map:
#             continue
#         for motif in motifs:
#             subseq = extract_subseq_from_seqs(motif, acc_seq_map[acc],
#                                               INIT_MOTIF_LEN, FINAL_MOTIF_LEN)
#             if subseq:
#                 assert len(subseq) == FINAL_MOTIF_LEN
#                 aligned_sequences.append(subseq)
#
#     matrix_counter = np.zeros([FINAL_MOTIF_LEN, 20], dtype=int)
#     for subseq in aligned_sequences:
#         if matrix_counter is None:
#         for i, char in enumerate(subseq):
#             try:
#                 AA_index = generic.AA1_TO_I[char]
#             except KeyError:
#                 # Key not found for some reason
#                 continue
#             matrix_counter[i, AA_index] += 1
#
#     build_meme_from_aligned._write_matrix_file(matrix_counter, matrix_file)
#     generic.quit_if_missing(matrix_file)
#     generic.quit_if_missing(composition)
#     command = f"{paths.MATRIX_2_MEME_EXEC} -protein -bg" \
#               f" {composition} < {matrix_file} > {output}"
#     subprocess.run(command, shell=True)
#     generic.quit_if_missing(output)
#     # end
#     # to check logo of this, run ./meme/ceqlogo -i1 output_meme_mine.txt -o
#     # logo.eps -f EPS
#     # output as matrix, and number of Kmatches.
#     # Number of Kmatches is approximately number of aligned seqs
#     return
#
#
# def main():
#     """
#     Two possible workflows at this point:
#     1. Prosite
#         1. Download aligned_seqs
#         2. Run run_prosite_aligned_cropped, output as number of Kmatches,
#         and input matrix for converge.
#         3. Encode into output_matrix_binary using the output_matrix.
#         4. Run converge using that output_matrix_binary as seed_alignment, on uniprot.
#         Manually update Kmatches in converge code.
#         5. With converge output and composition.txt, run conv_to_meme to
#         convert to meme.
#         6. Run ceqlogo on meme to get profile.
#
#     2. NBDB:
#         1. Download NBDB matrix file
#         2. Manually get number of Kmatches, from length of matching seqs.
#         3. Convert this NBDB matrix file to converge seed matrix file,
#         using convert_nbdb_matrix_to_conv_encoder()
#         4. See 3-6 above.
#
#     3. todo: remember that igor also wants one that's wider, at length 50.
#     One for flat pssm for the outer region, and one for included region.
#     """
#     # todo: build a composition file based on uniprot, then just use that.
#
#
#
#
#     # For GxGGxG: 2602
#     # For GxGxxG: 4077
#     # For GxxGxG: 4872
#     # todo: check how dependent Kmatches here is. Also, since pdb with
#     #  different cid but same pdb_id are probably the same thing, should the
#     #  Kmatches then be based on number of unique pdb_id rather than all...?
#     # logs.set_logging_level()
#     # matrix = convert_nbdb_matrix_to_conv_encoder(os.path.join(
#     #     paths.ROOT, "GxGxxG_pssm.txt"), 4077)
#     # with open("output_matrix.txt", "w") as file:
#     #     for line in matrix:
#     #         for value in line:
#     #             file.write(str(value) + " ")
#     #         file.write("\n")
#     from meme_suite import meme_interface
#     seq_file = os.path.join(paths.ROOT, 'efhand_seq.fasta')
#     motif_len = 15
#     meme_folder = paths.MEME_MAST_FOLDER
#     num_p = 7
#     from utils.clean_fasta_alphabet import screen
#     screen(fasta_path=seq_file, output_path=seq_file)
#     import filter_seqs
#     filter_seqs.delete_short_seqs(seq_file, threshold=15)
#     meme_interface.run_meme(seq_file, motif_len, meme_folder, num_p)
#     meme_txt = os.path.join(meme_folder, "meme.txt")
#
#     meme_interface.run_mast(meme_txt, seq_file, meme_folder)
#     mast_txt_path = os.path.join(meme_folder, 'mast.txt')
#     acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
#                                                                motif_len)
#
#     acc_seq_map = get_pname_seq.parse_raw(seq_file)
#     aligned_sequences = []
#     for acc, motifs in acc_motif_map.items():
#         if acc not in acc_seq_map:
#             continue
#         seq = acc_seq_map[acc]
#         if motif_len > 49:
#             print("motif_len bigger than 50. Exiting.")
#             raise Exception
#         left_spacing = 25 - motif_len // 2
#         right_spacing = 25 + motif_len // 2
#         for motif in motifs:
#             if len(seq) < motif + right_spacing:
#                 continue
#             if motif < left_spacing:
#                 continue
#             aligned_seq = seq[motif - left_spacing:motif + right_spacing]
#             assert len(aligned_seq) == 50
#             aligned_sequences.append(aligned_seq)
#
#     # if using meme and we need meme.txt, then run this:
#     # with open("aligned_seq_mine.txt", "w") as file:
#     #     for seq in aligned_sequences:
#     #         file.write(">RAND\n")
#     #         file.write(seq + "\n")
#     # composition_file = os.path.join(paths.ROOT, "uniprot_full_composition.txt")
#     # build_composition.build("uniprot_full.fasta", composition_file)
#     # build_meme_from_aligned.build("aligned_seq_mine.txt",
#     #                               "output_meme_mine.txt", composition_file)
#     # end
#
#     # to check logo of this, run ./meme/ceqlogo -i1 output_meme_mine.txt -o
#     # logo.eps -f EPS
#     # output as matrix, and number of Kmatches.
#     # Number of Kmatches is approximately number of aligned seqs
#     file_folder = os.path.join(paths.ROOT, "Reports", "19082019", "PDOC00706")
#
#     alphabets = dict()
#     AA_values = sorted(set(generic.AA3_to_AA1.values()))
#     assert len(AA_values) == 20
#     for i, AA in enumerate(AA_values):
#         alphabets[AA] = i
#     output_matrix = np.zeros((50, 20), dtype=int)
#     for seq in aligned_sequences:
#         for i, char in enumerate(seq):
#             if char in alphabets:
#                 output_matrix[i][alphabets[char]] += 1
#     # print(output_matrix)
#     output_matrix_path = os.path.join(file_folder, "output_matrix.txt")
#
#     with open(output_matrix_path, 'w') as file:
#         for line in output_matrix:
#             for value in line:
#                 file.write(str(value) + " ")
#             file.write("\n")
#     print(f"Kmatches is {len(aligned_sequences)}.")
#     return
#
#
#     sys.exit()
#     file_folder = os.path.join(paths.ROOT, "Reports", "19082019", "PDOC00706")
#
#     full_seq_path = os.path.join(file_folder, "full_seq.fasta")
#     aligned_seq_path = os.path.join(file_folder, "aligner.txt")
#     output_matrix_path = os.path.join(file_folder, "output_matrix.txt")
#     # preprocess.run_prosite_aligned_cropped_flat(full_seq_path, aligned_seq_path,
#     #                                        output_matrix_path)
#     # preprocess.run_prosite_aligned_cropped(full_seq_path, aligned_seq_path,
#     #                                             output_matrix_path)
#     # import sys
#     # sys.exit()
#
#     # 1. move output_matrix.txt to converge_encoder. Run converge_encoder.
#     # 2. Move input_matrix_binary into conv_rewrite. Take note of Kmatches
#     # here, and edit in single_Kmatches in conv_rewrite. Run conv_rewrite.
#     # 3. Move output.1.matrix and composition.txt back to Reports folder,
#     # run the rest below.
#
#
#
#
#     # After that
#     from converge.conv_interface import convert_conv_to_meme_full_num
#     converge_output_path = os.path.join(file_folder, "output.1.matrix")
#     composition_path = os.path.join(file_folder, "composition.txt")
#     output_meme_path = os.path.join(file_folder, "meme.txt")
#     convert_conv_to_meme_full_num(converge_output_path, composition_path,
#                                   output_meme_path)
#
#     logo_path = os.path.join(file_folder, "logo.eps")
#     import subprocess
#     ceqlogo_exec = os.path.join(paths.ROOT, "src", "meme_suite", "meme",
#                                 "ceqlogo")
#     command = f"{ceqlogo_exec} -i1 {output_meme_path} -o {logo_path} -f EPS"
#     subprocess.run(command, shell=True)
#
#     tmp_composition_path = os.path.join(file_folder, "compos_align.txt")
#     aligned_meme_path = os.path.join(file_folder, "aligned_meme.txt")
#     build_composition.build(aligned_seq_path, tmp_composition_path)
#     build_meme_from_aligned.build(aligned_seq_path, aligned_meme_path,
#                                   tmp_composition_path)
#     aligned_logo_path = os.path.join(file_folder, "aligned_logo.eps")
#     command = f"{ceqlogo_exec} -i1 {aligned_meme_path} -o {aligned_logo_path} -f EPS"
#     subprocess.run(command, shell=True)
#
#
#     # # with open(os.path.join(paths.ROOT, "motif_pos.pkl"), 'rb') as file:
#     # #     print(len(pickle.load(file)))
#     # preprocess.run_prosite_aligned_cropped(paths.PROSITE_ENOLASE_SEQS,
#     #                                paths.PROSITE_ALIGNED_SEQS,
#     #                                os.path.join(paths.ROOT,
#     #                                             "output_motif_pos_new.pkl"))
#
#
#
#     preprocess.run_prosite_mast(paths.PROSITE_EXTRACT, 30,
#                                 paths.REF_MEME_TXT,
#                      os.path.join(paths.ROOT, "prosite_mast_motif_pos.pkl"))
#     preprocess.run_prosite_aligned(paths.PROSITE_ENOLASE_SEQS,
#                                    paths.PROSITE_ALIGNED_SEQS,
#                                    os.path.join(paths.ROOT,
#                                                 "output_motif_pos.pkl"))
#
#     # seq_path = os.path.join(paths.USER_INPUT, "aligner.txt")
#     # matrix_ordered = matrix_builder(seq_path)
#     # write_matrix_file(matrix_ordered,
#     #                   os.path.join(paths.ROOT, "mine__.txt"))
#     # print(matrix_ordered)
#     # import sys
#     # sys.exit()
#     # logs.set_logging_level()
#     #
#     # from utils import get_pname_seq
#     #
#     # # pname_seq_map = get_pname_seq.parse(seq_path)
#     # from meme_suite import meme_interface
#     # meme_interface.create_meme_from_aligned(seq_path, 14, paths.MEME_MAST_FOLDER, num_p=7)
#     # meme_txt = os.path.join(paths.MEME_MAST_FOLDER, "meme.txt")
#
#
#     # logs.set_logging_level()
#     # extract_path = paths.PROSITE_EXTRACT
#     # motif_len = 13
#     # output = paths.PID_PDB_MAP
#     # num_p = 7
#     # preprocess.run_prosite_meme(extract_path, motif_len, output, num_p)
#
#     # extract_path = paths.PROSITE_EXTRACT
#     # motif_len = 13
#     # ref_meme_txt = paths.REF_MEME_TXT
#     # print(ref_meme_txt)
#     # output = paths.PID_PDB_MAP
#     # preprocess.run_prosite_mast(extract_path, motif_len, ref_meme_txt, output)
#     # with open(output, 'rb') as file:
#     #     print(pickle.load(file))
#     #
#     # extract_path = paths.IONCOM_EXTRACT
#     # motif_len = 13
#     # ref_meme_txt = paths.REF_MEME_TXT
#     # output = paths.PID_PDB_MAP
#     # preprocess.run_ioncom_mast(extract_path, motif_len, ref_meme_txt, output)
#     pass


# if __name__ == "__main__":
#     main()
seq_aligned_converge()
"""
Use this to extract meme matrix and merge them. 

def extract_meme_matrix(file):
    meme_matrix = []
    start = False
    width = None
    nsites = None
    with open(file, "r") as file:
        for line in file:
            if line.startswith("letter-probability"):
                start = True
                width = int(re.search("w= ([0-9]+)", line).group(1))
                nsites = int(re.search("nsites= ([0-9]+)", line).group(1))
                continue
            if not line.strip():
                continue
            if start:
                meme_matrix.append(list(map(float, line.strip().split(" "))))
    meme_matrix = np.array(meme_matrix)
    assert meme_matrix.shape[0] == width, meme_matrix.shape
    return (nsites, meme_matrix)
from config import paths

def meme_merger(file1, file2, offset_1_by=None, offset_2_by=None):
    nsites1, meme_matrix1 = extract_meme_matrix(file1)
    nsites2, meme_matrix2 = extract_meme_matrix(file2)
    meme_matrix2 /= 2
    meme_matrix1 /= 2
    if offset_1_by:
        if offset_1_by > 0:
            meme_matrix2[:(-offset_1_by)] += meme_matrix1[offset_1_by:]
        else:
            meme_matrix2[abs(offset_1_by):] += meme_matrix1[:offset_1_by]
        for i, line in enumerate(meme_matrix2):
            meme_matrix2[i] = line / sum(line)
        output_matrix = meme_matrix2

    elif offset_2_by:
        if offset_2_by > 0:
            meme_matrix1[:(-offset_2_by)] += meme_matrix2[offset_2_by:]
        else:
            meme_matrix1[abs(offset_2_by):] += meme_matrix2[:offset_2_by]
        for i, line in enumerate(meme_matrix1):
            meme_matrix1[i] = line / sum(line)
        output_matrix = meme_matrix1
    else:
        meme_matrix2 += meme_matrix1
        output_matrix = meme_matrix2
    return min(nsites1, nsites2), output_matrix


nsite, output_matrix = meme_merger(os.path.join(paths.ROOT,
                                           "output_meme_gxgxxg.txt"),
            os.path.join(paths.ROOT, "output_meme_gxxgxg.txt"), offset_2_by=1)


def write_matrix(output_matrix, nsite, output_path):
    with open(output_path, 'w') as file:
        file.write("MEME version 4\n\n")
        file.write("ALPHABET= ACDEFGHIKLMNPQRSTVWY\n\n")
        file.write("Background letter frequencies\n")
        file.write("A 0.0826 C 0.0138 D 0.0546 E 0.0673 F 0.0387 G 0.0708 H "
                   "0.0228 I 0.0592 K 0.0582 L 0.0966\n")
        file.write("M 0.0242 N 0.0406 P 0.0473 Q 0.0393 R 0.0554 S 0.0663 T "
                   "0.0536 V 0.0687 W 0.011\n")
        file.write("Y 0.0292\n\n")
        file.write("MOTIF MEME-1\n")
        file.write(f"letter-probability matrix: alength= 20 w="
                   f" {output_matrix.shape[0]} nsites= {nsite} "
                   f"E= 0.000\n")
        for line in output_matrix:
            for item in line:
                file.write("{:6f} ".format(item))
            file.write("\n")


output_path = os.path.join(paths.ROOT, 'written_meme.txt')
write_matrix(output_matrix, nsite, output_path)

"""









