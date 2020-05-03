import sys
import find_cid_from_pname
from utils import logs

from config import paths
import preprocess
import os
import re
from collections import Counter, OrderedDict
import time

import traceback
import logging

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



def seq_aligned_converge_testing():
    '''
    16/10/2019
    Testing to see if my conv_rewrite gives same output as the original
    converge. Mainly in tracking of maxS.

    To do this, we'll take aligned_sequences, run make_meme on it, then run
    meme_to_conv, then use that conv as seed_matrix, and run conv_orig on the
    full uniprot data.
    '''
    # aligned sequences need to have no gap
    # todo: test for that?

    # Known matched sequences
    SEQ_FILE = os.path.join(paths.ROOT, 'mg_seqs.fasta')
    ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'mg_aligned.txt')
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
            for value in line[:-1]:
                file.write(str(value) + ",")
            file.write(str(line[-1]) + "\n")
    print(f"Kmatches is {len(aligned_sequences)}.")

    #     next, we run converge.
    #     assumption that proteome, blosum is already encoded
    from converge import conv_interface

    OUTPUT_MEME = os.path.join(paths.ROOT, "output_meme.txt")

    conv_interface.run(OUTPUT_MATRIX, FINAL_MOTIF_LEN, len(aligned_sequences),
                       OUTPUT_MEME)


def main():
    pdb_seq_file = paths.RCSB_SEQS_FASTA
    clean_fasta_alphabet.screen(pdb_seq_file)

    uniprot_full = os.path.join(paths.USER_INPUT, "uniprot_full.fasta")
    composition = os.path.join(paths.INTERNAL, "uniprot_full_composition.txt")
    build_composition.build(uniprot_full, composition)

    print("Start mg_dxdxxd")
    timestamp = time.time()
    mg_seq = os.path.join(paths.ROOT, 'mg_seqs.fasta')
    mg_align = os.path.join(paths.ROOT, 'mg_aligned.txt')
    mg_output_matrix = os.path.join(paths.ROOT,
                                           'mg_output_matrix.txt')
    mg_numseq = build_output_matrix_from_aligned(mg_seq, mg_align,
                                     mg_output_matrix, composition)
    mg_motif_pos = os.path.join(paths.ROOT, "mg_motif_pos.txt")
    get_motif_pos_from_output_matrix(mg_output_matrix, mg_numseq,
                                     pdb_seq_file, mg_motif_pos)
    print("End mg_dxdxxd: {}".format(time.time() - timestamp))

    print("Start efhand")
    timestamp = time.time()
    efhand_seqs = os.path.join(paths.ROOT, 'efhand_aligned_seqs.fasta')
    efhand_align = os.path.join(paths.ROOT, 'efhand_aligned.txt')
    efhand_output_matrix = os.path.join(paths.ROOT, 'efhand_output_matrix.txt')
    efhand_numseq = build_output_matrix_from_aligned(efhand_seqs, efhand_align,
                                     efhand_output_matrix, composition)
    efhand_motif_pos = os.path.join(paths.ROOT, "efhand_motif_pos.txt")
    get_motif_pos_from_output_matrix(efhand_output_matrix, efhand_numseq,
                                     pdb_seq_file, efhand_motif_pos)
    print("End efhand: {}".format(time.time() - timestamp))

    # num_seqs don't really matter here because of normalisation
    # we cannot use numseq from NBDB file, cos it's too big. So just use a
    # random number, it doesn't really matter.
    print("Start GxGGxG")
    timestamp = time.time()
    GxGGxG_pssm = os.path.join(paths.ROOT, "GxGGxG_pssm.txt")
    GxGGxG_numseq = 2000
    GxGGxG_output_matrix = os.path.join(paths.ROOT, 'GxGGxG_output_matrix.txt')
    build_output_matrix_from_pssm(GxGGxG_pssm, GxGGxG_numseq, GxGGxG_output_matrix)
    GxGGxG_motif_pos = os.path.join(paths.ROOT, "GxGGxG_motif_pos.txt")
    get_motif_pos_from_output_matrix(GxGGxG_output_matrix, GxGGxG_numseq,
                                     pdb_seq_file, GxGGxG_motif_pos)
    print("End GxGGxG: {}".format(time.time() - timestamp))

    print("Start GxGxxG")
    GxGxxG_pssm = os.path.join(paths.ROOT, "GxGxxG_pssm.txt")
    GxGxxG_numseq = 2000
    GxGxxG_output_matrix = os.path.join(paths.ROOT, 'GxGxxG_output_matrix.txt')
    build_output_matrix_from_pssm(GxGxxG_pssm, GxGxxG_numseq, GxGxxG_output_matrix)
    GxGxxG_motif_pos = os.path.join(paths.ROOT, "GxGxxG_motif_pos.txt")
    get_motif_pos_from_output_matrix(GxGxxG_output_matrix, GxGxxG_numseq,
                                     pdb_seq_file, GxGxxG_motif_pos)
    print("End GxGxxG: {}".format(time.time() - timestamp))

    print("Start GxxGxG")
    timestamp = time.time()
    GxxGxG_pssm = os.path.join(paths.ROOT, "GxxGxG_pssm.txt")
    GxxGxG_numseq = 2000
    GxxGxG_output_matrix = os.path.join(paths.ROOT, 'GxxGxG_output_matrix.txt')
    build_output_matrix_from_pssm(GxxGxG_pssm, GxxGxG_numseq, GxxGxG_output_matrix)
    GxxGxG_motif_pos = os.path.join(paths.ROOT, "GxxGxG_motif_pos.txt")
    get_motif_pos_from_output_matrix(GxxGxG_output_matrix, GxxGxG_numseq,
                                     pdb_seq_file, GxxGxG_motif_pos)
    print("End GxxGxG: {}".format(time.time() - timestamp))


def build_output_matrix_from_aligned(selected_seqs, aligned_seqs, output,
                                     composition=None):
    """
    # SEQ_FILE = os.path.join(paths.ROOT, 'mg_seqs.fasta')
    # ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'mg_aligned.txt')

    # SEQ_FILE = os.path.join(paths.ROOT, 'efhand_aligned_seqs.fasta')
    # ALIGNED_SEQ_FILE = os.path.join(paths.ROOT, 'efhand_aligned.txt')
    :return:
    """
    # Length of actual motif
    FINAL_MOTIF_LEN = 30
    # composition_file should be provided, can be generated as below:
    if composition is None:
        composition = os.path.join(paths.ROOT, "uniprot_full_composition.txt")
        build_composition.build(os.path.join(paths.ROOT, "uniprot_full.fasta"), \
                                composition)

    # Intermediate paths
    MEME_FOLDER = paths.MEME_MAST_FOLDER
    MEME_TXT = paths.TMP_FILE_TEMPLATE.format(1)
    CLEANED_SEQ = paths.TMP_FILE_TEMPLATE.format(2)

    INIT_MOTIF_LEN = preprocess._get_motif_len_from_aligned(aligned_seqs)
    clean_fasta_alphabet.screen(fasta_path=selected_seqs, output_path=CLEANED_SEQ)
    filter_seqs.delete_short_seqs(CLEANED_SEQ, threshold=FINAL_MOTIF_LEN)

    build_meme_from_aligned.build(aligned_seqs, MEME_TXT, composition)

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

    output_matrix = np.zeros((FINAL_MOTIF_LEN, 20), dtype=int)
    for seq in aligned_sequences:
        for i, char in enumerate(seq):
            if char in generic.AA1_TO_I:
                output_matrix[i][generic.AA1_TO_I[char]] += 1

    with open(output, 'w') as file:
        for line in output_matrix:
            for value in line:
                file.write(str(value) + " ")
            file.write("\n")
    return len(aligned_sequences)


def build_output_matrix_from_pssm(pssm_file, num_seqs, output):
    """
    For nbdb pssm files
    GxGGxG_pssm: 4000
    """
    # composition_file should be provided, can be generated as below:
    output_matrix = convert_nbdb_matrix_to_conv_encoder(pssm_file, num_seqs)
    with open(output, 'w') as file:
        for line in output_matrix:
            for value in line:
                file.write(str(value) + " ")
            file.write("\n")


def get_motif_pos_from_output_matrix(output_matrix, num_seqs,
                                     pdb_seq_file, output, motif_len=30):
    from converge import conv_interface

    OUTPUT_MEME = os.path.join(paths.ROOT, "output_meme.txt")

    conv_interface.run(output_matrix, motif_len, num_seqs, OUTPUT_MEME)

    meme_interface.run_mast(OUTPUT_MEME, pdb_seq_file, paths.MEME_MAST_FOLDER)
    mast_txt_path = os.path.join(paths.MEME_MAST_FOLDER, 'mast.txt')

    pdb_cid_motif_raw = _get_motif_diagram_mast_pdb(mast_txt_path)

    pdb_cid_seq = dict()
    for i, (pdb_id, cid) in enumerate(pdb_cid_motif_raw.keys()):
        if not i % 10:
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
    with open(paths.TMP_FILE_TEMPLATE.format("tmp_pdb_cid_seq.pkl"), 'wb') as \
            file:
        pickle.dump(pdb_cid_seq, file, -1)
    pdb_cid_seq = OrderedDict(sorted(pdb_cid_seq.items()))
    pdb_seq_direct_from_pdb_files = paths.TMP_FILE_TEMPLATE.format(12)
    with open(pdb_seq_direct_from_pdb_files, 'w') as file:
        for (pdb_id, cid), seq in pdb_cid_seq.items():
            file.write(f">{pdb_id}_{cid}\n")
            file.write(seq + "\n")
    clean_fasta_alphabet.screen(pdb_seq_file)
    filter_seqs.delete_short_seqs(pdb_seq_direct_from_pdb_files, motif_len)
    meme_interface.run_mast(OUTPUT_MEME, pdb_seq_direct_from_pdb_files,
                            paths.MEME_MAST_FOLDER)
    mast_txt_path = os.path.join(paths.MEME_MAST_FOLDER, 'mast.txt')

    pdb_cid_motif_raw = _get_motif_diagram_mast_pdb(mast_txt_path)
    pdb_cid_motif_map = meme_interface._adjust_motif_diagram(pdb_cid_motif_raw,
                                                             motif_len)

    motif_positions = defaultdict(dict)
    for (pdb_id, cid), motif_pos in pdb_cid_motif_map.items():
        motif_positions[pdb_id]['sno_markers'] = motif_pos
        motif_positions[pdb_id]['cid'] = cid
    motif_positions = OrderedDict(sorted(motif_positions.items()))

    with open(output, 'wb') as file:
        pickle.dump(motif_positions, file, -1)


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

if __name__ == "__main__":
    main()