# todo: take pssm from nbdb website, find out how many nsites that converge
#  to, convert into the necessary binaries, then run converge on uniprot
#  using that. Once that is obtained, same thing, run on ceqlogo to get
#  profile. Do the same for all 3 tasks.

# In order to do that though, we need to know the number of nsites used to
# create the pssms. So, need igor's data...?  Well, we can guess that it'll
# be 1300, or number of seqs in the .srh files. 

# This dataset has 560537 sequences, of which we start with a seed alignment
# of 912, and end with 1589 matches.

import logging
import os
import pickle
import shutil

import numpy as np

from config import paths
import extract_parser, download_pdb_files, create_seq_file, filter_seqs, \
    motif_finder, prosite_pdb_list, find_cid_from_pname
from meme_suite import meme_interface
from converge import conv_interface
from utils import generic, logs, get_pname_seq, build_composition, \
    build_meme_from_aligned, download_fasta_given_pdb, uniprot_id_converter
from collections import defaultdict

from pdb_component import pdb_interface

import sys
import re
from collections import OrderedDict
import random


def _get_motif_len_from_aligned(aligned_file):
    with open(aligned_file, 'r') as file:
        for line in file:
            if not line.startswith(">"):
                return len(line.strip())

def keep_only_acc(acc_list, seq_file, output):
    header_seq_map = generic.read_fasta(seq_file)
    acc_list = set(acc_list)
    headers = list(header_seq_map.keys())
    for header in headers:
        acc = header.split("|")[1]
        if acc not in acc_list:
            del header_seq_map[header]
    generic.write_fasta(header_seq_map, output)


def run_prosite_aligned_cropped(seq_file, aligned_seq_file, output,
                                storage_path=None):
    """
    15/8/2019
    First, we produce a meme.txt based on aligned seqs.
    Then, we run mast on this, on the list of known matching seqs, from prosite.
    Then, we extract a length of 30, and collate this into another meme.txt.
    Then, we convert this meme.txt into a table form, for conv_encoder to
    take over.
    Conv_encoder then converts it into a PSSM, and run converge on the whole
    uniref/uniprot data. # note: see if converge, when running on matrix,
    actually converges it or it just finds matching ones.
    Finally, we collate the converged profile, and the sequences that match
    to it. Build profile, present to Igor.
    :return:
    """

    generic.quit_if_missing(seq_file)
    generic.quit_if_missing(aligned_seq_file)
    generic.warn_if_exist(output)
    if storage_path is None:
        composition_file = paths.TMP_FILE_TEMPLATE.format('composition.txt')
        meme_txt = paths.TMP_FILE_TEMPLATE.format('meme_from_aligned.txt')
        meme_folder = paths.MEME_MAST_FOLDER
        cropped_seq_file = paths.TMP_FILE_TEMPLATE.format('cropped_seqs.fasta')
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        composition_file = os.path.join(storage_path, 'composition.txt')
        meme_txt = os.path.join(storage_path, 'meme_from_aligned.txt')
        meme_folder = os.path.join(storage_path, 'meme_mast_folder')
        cropped_seq_file = os.path.join(storage_path, 'cropped_seqs.fasta')
    if os.path.isfile(composition_file):
        os.remove(composition_file)
    build_composition.build(seq_file, composition_file)
    motif_len = _get_motif_len_from_aligned(aligned_seq_file)
    build_meme_from_aligned.build(aligned_seq_file, meme_txt, composition_file)
    meme_interface.run_mast(meme_txt, seq_file, meme_folder)
    mast_txt_path = os.path.join(meme_folder, 'mast.txt')
    acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
                                                               motif_len)

    acc_seq_map = get_pname_seq.parse_raw(seq_file)
    aligned_sequences = []
    for acc, motifs in acc_motif_map.items():
        if acc not in acc_seq_map:
            continue
        seq = acc_seq_map[acc]
        if motif_len > 49:
            print("motif_len bigger than 50. Exiting.")
            raise Exception
        left_spacing = 25 - motif_len//2
        right_spacing = 25 + motif_len // 2
        for motif in motifs:
            if len(seq) < motif + right_spacing:
                continue
            if motif < left_spacing:
                continue
            aligned_seq = seq[motif-left_spacing:motif+right_spacing]
            aligned_sequences.append(aligned_seq)
    with open("aligned_seq_mine.txt", "w") as file:
        for seq in aligned_sequences:
            file.write(">RAND\n")
            file.write(seq+"\n")
    build_meme_from_aligned.build("aligned_seq_mine.txt", "output_meme_mine.txt",
                                  composition_file)
    # to check logo of this, run ./meme/ceqlogo -i1 output_meme_mine.txt -o
    # logo.eps -f EPS
    # output as matrix, and number of Kmatches.
    # Number of Kmatches is approximately number of aligned seqs
    alphabets = dict()
    AA_values = sorted(set(generic.AA3_to_AA1.values()))
    assert len(AA_values) == 20
    for i, AA in enumerate(AA_values):
        alphabets[AA] = i
    output_matrix = np.zeros((50, 20), dtype=int)
    for seq in aligned_sequences:
        for i, char in enumerate(seq):
            if char in alphabets:
                output_matrix[i][alphabets[char]] += 1
    # print(output_matrix)
    with open(output, 'w') as file:
        for line in output_matrix:
            for value in line:
                file.write(str(value) + " ")
            file.write("\n")
    print(f"Kmatches is {len(aligned_sequences)}.")
    return



def run_prosite_aligned_cropped_flat(seq_file, aligned_seq_file, output,
                                storage_path=None):
    """
    15/8/2019
    Flatten the areas outside the aligned portions.
    :return:
    """

    generic.quit_if_missing(seq_file)
    generic.quit_if_missing(aligned_seq_file)
    generic.warn_if_exist(output)
    if storage_path is None:
        composition_file = paths.TMP_FILE_TEMPLATE.format('composition.txt')
        meme_txt = paths.TMP_FILE_TEMPLATE.format('meme_from_aligned.txt')
        meme_folder = paths.MEME_MAST_FOLDER
        cropped_seq_file = paths.TMP_FILE_TEMPLATE.format('cropped_seqs.fasta')
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        composition_file = os.path.join(storage_path, 'composition.txt')
        meme_txt = os.path.join(storage_path, 'meme_from_aligned.txt')
        meme_folder = os.path.join(storage_path, 'meme_mast_folder')
        cropped_seq_file = os.path.join(storage_path, 'cropped_seqs.fasta')
    if os.path.isfile(composition_file):
        os.remove(composition_file)
    build_composition.build(seq_file, composition_file)
    composition_alphabets = []
    composition_probs = []

    with open(composition_file, 'r') as file:
        for line in file:
            if re.match("[A-Z] ", line):
                composition_alphabets.append(line[0])
                composition_probs.append(float(line[2:]))

    build_meme_from_aligned.build(aligned_seq_file, meme_txt, composition_file)
    motif_len = _get_motif_len_from_aligned(aligned_seq_file)
    meme_interface.run_mast(meme_txt, seq_file, meme_folder)
    mast_txt_path = os.path.join(meme_folder, 'mast.txt')
    acc_motif_map = meme_interface.extract_motifs_mast_uniprot(mast_txt_path,
                                                               motif_len)

    acc_seq_map = get_pname_seq.parse_raw(seq_file)
    aligned_sequences = []
    for acc, motifs in acc_motif_map.items():
        if acc not in acc_seq_map:
            continue
        seq = acc_seq_map[acc]
        for motif in motifs:
            if len(seq) < motif + 25:
                continue
            if motif < 25:
                continue
            aligned_seq = seq[motif:motif+motif_len]
            # replace the others with composition at random
            padding = random.choices(composition_alphabets,
                                     weights=composition_probs, k=50-motif_len)
            aligned_seq = "".join(padding[:15]) + aligned_seq \
                          + "".join(padding[15:])
            assert len(aligned_seq) == 50
            aligned_sequences.append(aligned_seq)
    # with open("aligned_seq_mine.txt", "w") as file:
    #     for seq in aligned_sequences:
    #         file.write(">RAND\n")
    #         file.write(seq + "\n")
    # build_meme_from_aligned.build("aligned_seq_mine.txt",
    #                               "output_meme_mine.txt", composition_file)
    # to check logo of this, run ./meme/ceqlogo -i1 output_meme_mine.txt -o
    # logo.eps -f EPS
    # output as matrix, and number of Kmatches.
    # Number of Kmatches is approximately number of aligned seqs
    alphabets = dict()
    AA_values = sorted(set(generic.AA3_to_AA1.values()))
    assert len(AA_values) == 20
    for i, AA in enumerate(AA_values):
        alphabets[AA] = i
    output_matrix = np.zeros((50, 20), dtype=int)
    for seq in aligned_sequences:
        for i, char in enumerate(seq):
            if char in alphabets:
                output_matrix[i][alphabets[char]] += 1
    # print(output_matrix)
    with open(output, 'w') as file:
        for line in output_matrix:
            for value in line:
                file.write(str(value) + " ")
            file.write("\n")
    print(f"Kmatches is {len(aligned_sequences)}.")
    return


def run_prosite_aligned(seq_file, aligned_seq_file, output, storage_path=None):
    """
    We start off with a seq_file, aligned_seq_file.

    We derive matrix from aligned_seq_file
    We derive composition from seq_file
    We put both together into a meme.txt

    We take seq_file, and extract from it the relevant .pdb files and the
    corresponding cid. So acc+seq => pdb+cid

    We therefore select the seq files for which a corresponding .pdb+cid
    exists, and output it into a cropped_seqfile

    We run mast on this cropped_seqfile, to obtain motif_pos.

    Finally we output this motif_pos using the pdb+cid from before.

    Desired intermediate files:
    meme.txt from aligned_seq
    acc=>pdb+cid map
    cropped_seqfile.fasta
    acc=>motif_pos
    pdb+cid=>motif_pos
    """
    generic.quit_if_missing(seq_file)
    generic.quit_if_missing(aligned_seq_file)
    generic.warn_if_exist(output)
    if storage_path is None:
        composition_file = paths.TMP_FILE_TEMPLATE.format('composition.txt')
        meme_txt = paths.TMP_FILE_TEMPLATE.format('meme_from_aligned.txt')
        meme_folder = paths.MEME_MAST_FOLDER
        cropped_seq_file = paths.TMP_FILE_TEMPLATE.format('cropped_seqs.fasta')
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        composition_file = os.path.join(storage_path, 'composition.txt')
        meme_txt = os.path.join(storage_path, 'meme_from_aligned.txt')
        meme_folder = os.path.join(storage_path, 'meme_mast_folder')
        cropped_seq_file = os.path.join(storage_path, 'cropped_seqs.fasta')
    if os.path.isfile(composition_file):
        os.remove(composition_file)
    build_composition.build(seq_file, composition_file)
    build_meme_from_aligned.build(aligned_seq_file, meme_txt,
                                  composition_file)

    acc_seq_map = get_pname_seq.parse_raw(seq_file)
    acc_ids = list(acc_seq_map.keys())
    acc_pdb_map = uniprot_id_converter.convert("ACC", "PDB_ID", acc_ids)

    pdb_seq_map = dict()
    # because pdb => acc mapping may not be 1-1, we retain the original maps
    mapped_pdb_acc = dict()
    for acc_id, seq in acc_seq_map.items():
        if acc_id in acc_pdb_map:
            pdb_id = acc_pdb_map[acc_id]
            mapped_pdb_acc[pdb_id] = acc_id
            pdb_seq_map[pdb_id] = seq

    pdb_cid_map = find_cid_from_pname.find(pdb_seq_map)
    acc_pdb_cid_map = {mapped_pdb_acc[pdb]: (pdb, cid) for pdb, cid in
                           pdb_cid_map.items()}
    cropped_acc_list = list(acc_pdb_cid_map.keys())

    keep_only_acc(cropped_acc_list, seq_file, cropped_seq_file)
    meme_interface.run_mast(meme_txt, cropped_seq_file, meme_folder)
    mast_txt_path = os.path.join(meme_folder, 'mast.txt')
    acc_motif_map = meme_interface.extract_motifs_mast_uniprot(
        mast_txt_path, 14)
    acc_motif_map = motif_finder._delete_gapped_motifs_uniprot(acc_motif_map,
                                                               cropped_seq_file)
    pdb_motif_pos = defaultdict(dict)
    for acc, motif_pos in acc_motif_map.items():
        pdb_id, cid = acc_pdb_cid_map[acc]
        pdb_motif_pos[pdb_id]['sno_markers'] = motif_pos
        pdb_motif_pos[pdb_id]['cid'] = cid
    with open(output, 'wb') as file:
        pickle.dump(pdb_motif_pos, file, -1)
    if storage_path is None:
        shutil.move(composition_file, paths.TRASH)
        shutil.move(meme_txt, paths.TRASH)
        shutil.move(meme_folder, paths.TRASH)
        shutil.move(cropped_seq_file, paths.TRASH)
    return


def run_prosite_mast(extract_path, motif_len, ref_meme_txt, output,
                     pdb_folder=paths.PDB_FOLDER, storage_path=None):
    """
    :param extract_path: paths.IONCOM_EXTRACT
    :param motif_len: 13
    :param ref_meme_txt: paths.REF_MEME_TXT
    :param output: paths.PID_PDB_MAP
    :param pdb_folder:
    :param storage_path:
    :return:
    """
    generic.quit_if_missing(extract_path)
    generic.quit_if_missing(ref_meme_txt)
    generic.warn_if_exist(output)
    if storage_path is None:
        pname_cid_path = paths.PNAME_CID
        seq_path = paths.FULL_SEQS
        meme_folder = paths.MEME_MAST_FOLDER
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
        seq_path = os.path.join(storage_path, 'seqs.fasta')
        meme_folder = os.path.join(storage_path, 'meme_folder')

    parse_extract_prosite(extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, seq_path, pdb_folder)
    filter_seq_file(seq_path, threshold=31)
    find_motifs_mast(pname_cid_path, seq_path, ref_meme_txt, motif_len, output,
                     meme_folder)

    if storage_path is None:
        shutil.move(pname_cid_path, paths.TRASH)
        shutil.move(seq_path, paths.TRASH)
        shutil.move(meme_folder, paths.TRASH)


def run_prosite_meme(extract_path, motif_len, output, num_p=7,
                     pdb_folder=paths.PDB_FOLDER, storage_path=None):
    """
    :param extract_path: paths.PROSITE_EXTRACT
    :param motif_len: 13
    :param output: paths.PID_PDB_MAP
    """''
    generic.quit_if_missing(extract_path)
    generic.warn_if_exist(output)
    assert isinstance(num_p, int)
    assert num_p >= 1

    if storage_path is None:
        pname_cid_path = paths.PNAME_CID
        seq_path = paths.FULL_SEQS
        meme_folder = paths.MEME_MAST_FOLDER
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
        seq_path = os.path.join(storage_path, 'seqs.fasta')
        meme_folder = os.path.join(storage_path, 'meme_folder')

    parse_extract_prosite(extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, seq_path, pdb_folder)

    filter_seq_file(seq_path, threshold=31)

    find_motifs_meme(pname_cid_path, seq_path, motif_len, output,
                     meme_folder, num_p)
    if storage_path is None:
        shutil.move(pname_cid_path, paths.TRASH)
        shutil.move(seq_path, paths.TRASH)
        shutil.move(meme_folder, paths.TRASH)


def run_ioncom_mast(extract_path, motif_len, ref_meme_txt, output,
                    pdb_folder=paths.PDB_FOLDER, storage_path=None):
    """
    :param extract_path: paths.IONCOM_EXTRACT
    :param motif_len: 13
    :param ref_meme_txt: paths.REF_MEME_TXT
    :param output: paths.PID_PDB_MAP
    """
    generic.quit_if_missing(extract_path)
    generic.quit_if_missing(ref_meme_txt)
    generic.warn_if_exist(output)
    if storage_path is None:
        pname_cid_path = paths.PNAME_CID
        seq_path = paths.FULL_SEQS
        meme_folder = paths.MEME_MAST_FOLDER
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
        seq_path = os.path.join(storage_path, 'seqs.fasta')
        meme_folder = os.path.join(storage_path, 'meme_folder')

    parse_extract_ioncom(extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, seq_path, pdb_folder)

    filter_seq_file(seq_path, threshold=31)
    find_motifs_mast(pname_cid_path, seq_path, ref_meme_txt, motif_len, output,
                     meme_folder)
    if storage_path is None:
        shutil.move(pname_cid_path, paths.TRASH)
        shutil.move(seq_path, paths.TRASH)
        shutil.move(meme_folder, paths.TRASH)


def parse_extract_ioncom(input_file, pname_cid_path):
    """
    :param input_file: paths.IONCOM_EXTRACT
    :param pname_cid_path: paths.PNAME_CID
    """
    generic.quit_if_missing(input_file)
    pname_cid_map = extract_parser.parse_ioncom(input_file)
    _test_seq_cid_map(pname_cid_map)
    generic.warn_if_exist(pname_cid_path)
    with open(pname_cid_path, 'wb') as file:
        pickle.dump(pname_cid_map, file, -1)


def parse_extract_prosite(input_file, pname_cid_path):
    """
    :param input_file: paths.PROSITE_EXTRACT
    :param pname_cid_path: paths.PNAME_CID
    """
    generic.quit_if_missing(input_file)
    pname_cid_map = extract_parser.parse_prosite(input_file,
                                                 prosite_pdb_list.pdb_list)
    _test_seq_cid_map(pname_cid_map)
    generic.warn_if_exist(pname_cid_path)
    with open(pname_cid_path, 'wb') as file:
        pickle.dump(pname_cid_map, file, -1)


def download_pdb(pname_cid_path, pdb_folder=paths.PDB_FOLDER):
    """
    :param pname_cid_path: paths.PNAME_CID
    """
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    if not os.path.isdir(pdb_folder):
        logging.warning(f"PDB_folder in <{pdb_folder}> not found. Downloading "
                        f"from scratch takes a while.")
        os.mkdir(pdb_folder)
    download_pdb_files.download(pname_cid_map, pdb_folder)


def trim_pnames_based_on_pdb(pname_cid_path, pdb_folder=paths.PDB_FOLDER):
    """
    :param pname_cid_path: paths.PNAME_CID
    """
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    download_pdb_files.trim_pname_cid(pname_cid_map, pdb_folder)
    _test_seq_cid_map(pname_cid_map)
    with open(pname_cid_path, 'wb') as file:
        pickle.dump(pname_cid_map, file, -1)
    generic.quit_if_missing(pname_cid_path)


def create_seq(pname_cid_path, seq_path, pdb_folder=paths.PDB_FOLDER):
    """
    :param pname_cid_path: paths.PNAME_CID
    :param seq_path: paths.FULL_SEQS
    """
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    seqs = create_seq_file.extract_sequence(pdb_folder, pname_cid_map,
                                            AA3_to_AA1=generic.AA3_to_AA1)
    with open(seq_path, 'w') as file:
        file.writelines(seqs)
    create_seq_file.test_fasta_match_pdb(seq_path, pdb_folder, pname_cid_map,
                                         generic.AA3_to_AA1)


def filter_seq_file(seq_path, threshold=50):
    """
    :param seq_path: paths.FULL_SEQS
    """
    generic.quit_if_missing(seq_path)
    filter_seqs.delete_short_seqs(seq_path, threshold)


def find_motifs_meme(pname_cid_path, seq_file, motif_len, output,
                     meme_folder=paths.MEME_MAST_FOLDER, num_p=1):
    """
    :param pname_cid_path: paths.PNAME_CID
    :param seq_file: paths.FULL_SEQS
    :param motif_len: 13
    :param output: paths.MOTIF_POS
    """
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    motif_pos = motif_finder.find_meme(pname_cid_map,
                                       motif_len,
                                       num_p=num_p,
                                       meme_folder=meme_folder,
                                       seq_file=seq_file)
    generic.warn_if_exist(output)
    with open(output, 'wb') as file:
        pickle.dump(motif_pos, file, -1)


def find_motifs_mast_uniprot(pname_cid_path, seq_file, ref_meme_txt, motif_len,
                       output,
                     meme_folder=paths.MEME_MAST_FOLDER):
    """
    :param pname_cid_path: paths.PNAME_CID
    :param seq_file: paths.FULL_SEQS
    :param ref_meme_txt: paths.REF_MEME_TXT
    :param motif_len: 13
    :param output: paths.MOTIF_POS
    """
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    motif_pos = motif_finder.find_mast_uniprot(pname_cid_map, seq_file,
                                           ref_meme_txt,
                                       motif_len, meme_folder=meme_folder)
    generic.warn_if_exist(output)
    with open(output, 'wb') as file:
        pickle.dump(motif_pos, file, -1)


def find_motifs_mast(pname_cid_path, seq_file, ref_meme_txt, motif_len,
                     output, meme_folder=paths.MEME_MAST_FOLDER):
    """
    :param pname_cid_path: paths.PNAME_CID
    :param seq_file: paths.FULL_SEQS
    :param ref_meme_txt: paths.REF_MEME_TXT
    :param motif_len: 13
    :param output: paths.MOTIF_POS
    """
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    motif_pos = motif_finder.find_mast(pname_cid_map, seq_file,
                                       ref_meme_txt, motif_len,
                                       meme_folder=meme_folder)
    generic.warn_if_exist(output)
    with open(output, 'wb') as file:
        pickle.dump(motif_pos, file, -1)


def _test_seq_cid_map(seq_cid_map):
    for pname, cid in seq_cid_map.items():
        assert isinstance(pname, str)
        assert isinstance(cid, str)
        assert len(pname) == 4
        assert len(cid) == 1


"""

def run_converge(seq_path,
                 output,
                 binding_sites=None,
                 seed_seqs=None,
                 bash_exec=paths.BASH_EXEC,
                 num_p=7,
                 storage_path=None):
    "":param
seq_path: paths.FULL_SEQS:param
output: paths.MOTIF_POS:param
binding_sites: paths.IONCOM_BINDING_SITES:param
seed_seqs: paths.CONV_SEED_SEQS
""
assert isinstance(num_p, int)
assert num_p >= 1
assert os.path.isfile(seq_path)
if storage_path is None:
    pname_cid_path = paths.PNAME_CID
    conv_meme_file = paths.CONV_MEME_FILE
    meme_folder = paths.MEME_MAST_FOLDER
    motif_pos_path = paths.MOTIF_POS
    if not seed_seqs:
        seed_seqs = paths.CONV_SEED_SEQS
else:
    generic.quit_if_missing(storage_path, filetype='folder')
    pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
    conv_meme_file = os.path.join(storage_path, 'conv_meme.txt')
    meme_folder = os.path.join(storage_path, 'meme_folder')
    motif_pos_path = os.path.join(storage_path, 'motif_pos.pkl')
    if not seed_seqs:
        seed_seqs = os.path.join(storage_path, 'seed_seqs.fasta')

generic.warn_if_exist(conv_meme_file)
filter_seq_file(seq_path)

if binding_sites:
    generic.quit_if_missing(binding_sites)
    generic.warn_if_exist(seed_seqs)
    conv_interface.make_seed_seq(binding_sites, seed_seqs)

generic.quit_if_missing(seed_seqs)
conv_interface.run(seq_path, conv_meme_file, bash_exec, num_p, seed_seqs)

# meme_interface.run_mast(conv_meme_file, seq_path, meme_folder)
# mast_txt_path = os.path.join(meme_folder, 'mast.txt')
# motif_len = 30
# pname_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
#                                                      motif_len)
pname_seq_map = get_pname_seq.parse(seq_path)
download_pdb_files.trim_pname_cid(pname_seq_map, paths.PDB_FOLDER)
pname_cid_map = find_cid_from_pname.find(pname_seq_map, paths.PDB_FOLDER)
generic.warn_if_exist(pname_cid_path)
with open(pname_cid_path, 'wb') as file:
    pickle.dump(pname_cid_map, file, -1)

# Converge does not work now because it has multiple profiles in its
# meme.txt file, but we're only looking for one, in mast.

# find_motifs_mast(pname_cid_path, seq_path, conv_meme_file, 30,
#                  motif_pos_path, meme_folder=paths.MEME_MAST_FOLDER)
# with open(motif_pos_path, 'rb') as file:
#     print(pickle.load(file))
# print("\n")
# load_pdb_info(motif_pos_path, output)
# if storage_path is None:
#     shutil.move(conv_meme_file, paths.TRASH)
#     shutil.move(meme_folder, paths.TRASH)
#     shutil.move(motif_pos_path, paths.TRASH)
"""