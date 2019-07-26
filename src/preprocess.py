import logging
import os
import pickle
import shutil

from config import paths
import extract_parser, download_pdb_files, create_seq_file, filter_seqs, \
    motif_finder, prosite_pdb_list, loaders, find_cid_from_pname
from meme_suite import meme_interface
from converge import conv_interface
from utils import generic, logs, get_pname_seq


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
        motif_pos_path = paths.MOTIF_POS
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
        seq_path = os.path.join(storage_path, 'seqs.fasta')
        meme_folder = os.path.join(storage_path, 'meme_folder')
        motif_pos_path = os.path.join(storage_path, 'motif_pos.pkl')

    parse_extract_prosite(extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, seq_path, pdb_folder)

    filter_seq_file(seq_path, threshold=31)
    find_motifs_mast(pname_cid_path, seq_path, ref_meme_txt, motif_len,
                     motif_pos_path, meme_folder)

    load_pdb_info(motif_pos_path, output)
    if storage_path is None:
        shutil.move(pname_cid_path, paths.TRASH)
        shutil.move(seq_path, paths.TRASH)
        shutil.move(meme_folder, paths.TRASH)
        shutil.move(motif_pos_path, paths.TRASH)


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
        motif_pos_path = paths.MOTIF_POS
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
        seq_path = os.path.join(storage_path, 'seqs.fasta')
        meme_folder = os.path.join(storage_path, 'meme_folder')
        motif_pos_path = os.path.join(storage_path, 'motif_pos.pkl')

    parse_extract_prosite(extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, seq_path, pdb_folder)

    filter_seq_file(seq_path, threshold=31)

    find_motifs_meme(pname_cid_path, seq_path, motif_len, motif_pos_path,
                     meme_folder, num_p)
    load_pdb_info(motif_pos_path, output)
    if storage_path is None:
        shutil.move(pname_cid_path, paths.TRASH)
        shutil.move(seq_path, paths.TRASH)
        shutil.move(meme_folder, paths.TRASH)
        shutil.move(motif_pos_path, paths.TRASH)


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
        motif_pos_path = paths.MOTIF_POS
    else:
        generic.quit_if_missing(storage_path, filetype='folder')
        pname_cid_path = os.path.join(storage_path, 'pname_cid_map.pkl')
        seq_path = os.path.join(storage_path, 'seqs.fasta')
        meme_folder = os.path.join(storage_path, 'meme_folder')
        motif_pos_path = os.path.join(storage_path, 'motif_pos.pkl')

    parse_extract_ioncom(extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, seq_path, pdb_folder)

    filter_seq_file(seq_path, threshold=31)
    find_motifs_mast(pname_cid_path, seq_path, ref_meme_txt, motif_len,
                     motif_pos_path, meme_folder)
    load_pdb_info(motif_pos_path, output)
    if storage_path is None:
        shutil.move(pname_cid_path, paths.TRASH)
        shutil.move(seq_path, paths.TRASH)
        shutil.move(meme_folder, paths.TRASH)
        shutil.move(motif_pos_path, paths.TRASH)

def run_converge(seq_path,
                 output,
                 binding_sites=None,
                 seed_seqs=None,
                 bash_exec=paths.BASH_EXEC,
                 num_p=7,
                 storage_path=None):
    """
    :param seq_path: paths.FULL_SEQS
    :param output: paths.MOTIF_POS
    :param binding_sites: paths.IONCOM_BINDING_SITES
    :param seed_seqs: paths.CONV_SEED_SEQS
    """
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
    pname_cid_map = find_cid_from_pname.find(pname_seq_map)
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


def load_pdb_info(motif_pos_path, output):
    with open(motif_pos_path, 'rb') as file:
        motif_pos = pickle.load(file)
    pid_pdb_map = loaders.load_pdb_info(motif_pos)
    with open(output, 'wb') as file:
        pickle.dump(pid_pdb_map, file, -1)

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