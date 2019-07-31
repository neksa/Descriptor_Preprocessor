import logging
import os
import pickle
import shutil

import create_seq_file
from config import paths
import extract_parser
import filter_seqs
import motif_finder
from pdb_component import pdb_interface
import prosite_pdb_list
from utils import generic


def run_prosite_mast(extract_path, motif_len, ref_meme_txt, output,
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
    for pdb_id in pname_cid_map.keys():
        pdb_interface.download(pdb_id, silent=True)


def trim_pnames_based_on_pdb(pname_cid_path, pdb_folder=paths.PDB_FOLDER):
    """
    :param pname_cid_path: paths.PNAME_CID
    """
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    fnames = []
    for fname in os.listdir(pdb_folder):
        split_fname = fname.split(".")
        assert len(split_fname) == 2
        pname = split_fname[0]
        fnames.append(pname)
    fnames = set(fnames)
    pname_list = list(pname for pname in pname_cid_map.keys())
    for pname in pname_list:
        if pname not in fnames:
            del pname_cid_map[pname]
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


def filter_seq_file(seq_path, threshold=30):
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


def find_motifs_mast_uniprot(pname_cid_path, seq_file, ref_meme_txt,
                             motif_len, output,
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
                                               ref_meme_txt, motif_len,
                                               meme_folder=meme_folder)
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
