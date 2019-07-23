import datetime
import logging
import os
import pickle
import shutil

from config import paths
import extract_parser, download_pdb_files, create_seq_file, filter_seqs, \
    motif_finder, prosite_pdb_list
from meme_suite import meme_interface
from converge import conv_interface
from utils import generic, logs

# Using Converge
#############################################################
# First, we have to derive the converge seed sequences somehow. This should
# be stored in paths.CONV_SEED_SEQS.

# Possible ways are:
# 1. Manually
# 2. By random sampling of a sequence file
# 3. Based on Ioncom binding site

# Next, with a provided seqs_file, and seed_seqs_file, we run converge.

# Given a set of pnames, we need to
# 1. download the relevant .pdb files
# 2. if not inside, say so. If inside, give the cid. 

# Then, converge's output is converted to meme_format.

# Then, we run mast on that meme file, on the same or a different set of
# sequences

# Finally, we identify the motif positions based on the mast output,
# and pickle.dump that.


def run_all(process='meme', source='prosite', motif_len=13, num_p=7,
            extract_path=None, pdb_folder=paths.PDB_FOLDER,
            ref_meme_txt=paths.REF_MEME_TXT, output=paths.MOTIF_POS,
            storage_path=None):
    assert isinstance(num_p, int)
    assert num_p >= 1
    assert process in ('meme', 'mast')
    assert source in ('prosite', 'ioncom')

    if extract_path is None:
        if source == 'prosite':
            extract_path = paths.PROSITE_EXTRACT
        else:
            extract_path = paths.IONCOM_EXTRACT

    if storage_path is None:
        pname_cid_path = paths.PNAME_CID
        seq_path = paths.FULL_SEQS
        meme_folder = paths.MEME_MAST_FOLDER

    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
        seq_path = os.path.join(storage_path, 'seqs.fasta')
        meme_folder = os.path.join(storage_path, 'meme_folder')

    parse_extracts(source, extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, seq_path, pdb_folder)

    filter_seq_file(seq_path, threshold=31)
    find_motifs(process, motif_len, pname_cid_path, ref_meme_txt,
                meme_folder, seq_path, output, num_p)


# def run_converge(seq_path=paths.FULL_SEQS, output_path=paths.MOTIF_POS,
#                  binding_sites=paths.IONCOM_BINDING_SITES,
#                  seed_seqs_path=paths.CONV_SEED_SEQS,
#                  bash_exec=paths.BASH_EXEC,
#                  num_p=7, conv_meme_file=paths.CONV_MEME_FILE,
#                  mast_meme_folder=paths.MEME_MAST_FOLDER,
#                  delete_intermediate=False, storage_path=None):
#     assert isinstance(num_p, int)
#     assert num_p >= 1
#     assert os.path.isfile(seq_path)
#     generic.warn_if_exist(conv_meme_file)
#     filter_seq_file(seq_path)
#     conv_interface.run(seq_path, conv_meme_file, binding_sites, bash_exec,
#                        num_p, seed_seqs_path)
#
#     meme_interface.run_mast(conv_meme_file, seq_path, mast_meme_folder)
#     mast_txt_path = os.path.join(mast_meme_folder, 'mast.txt')
#     motif_len = 30
#     seq_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
#     motif_len)
#
#     generic.warn_if_exist(output_path)
#     with open(output_path, 'wb') as file:
#         pickle.dump(seq_motif_map, file, -1)
#     assert os.path.isfile(output_path)
#     if delete_intermediate:
#         os.remove(conv_meme_file)
#         shutil.rmtree(mast_meme_folder)
#     elif storage_path:
#         generic.warn_if_exist(storage_path, filetype='folder')
#         if not os.path.isdir(storage_path):
#             os.mkdir(storage_path)
#         shutil.move(conv_meme_file, storage_path)
#         shutil.move(mast_meme_folder, storage_path)
#     return


def parse_extracts(source, input_file_path, pname_cid_path):
    """
    :param source: 'prosite'
    :param input_file_path: extract_path = paths.IONCOM_EXTRACT
    :param pname_cid_path: paths.PNAME_CID
    :return: pname_cid_path
    """
    assert source in ('prosite', 'ioncom')
    if input_file_path == None:
        if source == 'prosite':
            input_file_path = paths.PROSITE_EXTRACT
        else:
            input_file_path = paths.IONCOM_EXTRACT
    generic.quit_if_missing(input_file_path)
    if source == 'prosite':
        pname_cid_map = extract_parser.parse_prosite(input_file_path,
                                                     prosite_pdb_list.pdb_list)
    else:
        pname_cid_map = extract_parser.parse_ioncom(input_file_path)
    _test_seq_cid_map(pname_cid_map)
    generic.warn_if_exist(pname_cid_path)
    with open(pname_cid_path, 'wb') as file:
        pickle.dump(pname_cid_map, file, -1)


def download_pdb(pname_cid_path, pdb_folder=paths.PDB_FOLDER):
    """
    :param pname_cid_path: paths.PNAME_CID
    :param pdb_folder: paths.PDB_FOLDER
    :return:
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
    :param pdb_folder: paths.PDB_FOLDER
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
    :param pdb_folder:
    :return: seq_path
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
    :param threshold: 30
    :return: seq_path
    """
    generic.quit_if_missing(seq_path)
    filter_seqs.delete_short_seqs(seq_path, threshold)


def find_motifs(process, motif_len, pname_cid_path,
                ref_meme_txt=paths.REF_MEME_TXT,
                meme_folder=paths.MEME_MAST_FOLDER,
                seq_file=paths.FULL_SEQS, output=paths.MOTIF_POS, num_p=1):
    """
    :param process: 'meme'
    :param motif_len: 13
    :param pname_cid_path: paths.PNAME_CID
    :param ref_meme_txt:
    :param mast_meme_folder:
    :param seq_file:
    :param output:
    :param num_p:
    :return:
    """
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    assert process in ('mast', 'meme')
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    motif_pos = motif_finder.find(pname_cid_map, motif_len, process=process,
                                  num_p=num_p, ref_meme_txt=ref_meme_txt,
                                  mast_meme_folder=mast_meme_folder,
                                  seq_file=seq_file)
    generic.warn_if_exist(output)
    with open(output, 'wb') as file:
        pickle.dump(motif_pos, file, -1)


def _test_seq_cid_map(seq_cid_map):
    for pname, cid in seq_cid_map.items():
        assert isinstance(pname, str)
        assert isinstance(cid, str)
        assert len(pname) == 4
        assert len(cid) == 1
    return True