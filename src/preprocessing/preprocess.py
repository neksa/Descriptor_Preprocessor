import datetime
import logging
import os
import pickle
import shutil

from config import paths
from preprocessing import extract_parser, download_pdb_files, \
    create_seq_file, filter_seqs, motif_finder, \
    prosite_pdb_list
from preprocessing.converge import conv_interface
from utils import generic, logs, run_mast_on_meme

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


def run_converge(seq_path=paths.FULL_SEQS,
                 output_path=paths.MOTIF_POS,
                 binding_sites=paths.IONCOM_BINDING_SITES,
                 seed_seqs_path=paths.CONV_SEED_SEQS,
                 bash_exec=paths.BASH_EXEC,
                 num_p=7,
                 conv_meme_file=paths.CONV_MEME_FILE,
                 mast_meme_folder=paths.MEME_MAST_FOLDER,
                 delete_intermediate=False,
                 storage_path=None):
    assert isinstance(num_p, int)
    assert num_p >= 1
    assert os.path.isfile(seq_path)
    generic.warn_if_exist(conv_meme_file)
    filter_seq_file(seq_path)
    conv_interface.run(seq_path, conv_meme_file, bash_exec, num_p,
                       seed_seqs_path)
    run_mast_on_meme.create(conv_meme_file, mast_meme_folder, seq_path)
    seq_motif_map = motif_finder._build_seq_motif_map("mast", mast_meme_folder,
                                                      seq_path,
                                      num_p=1, ref_meme_txt=conv_meme_file)
    generic.warn_if_exist(output_path)
    with open(output_path, 'wb') as file:
        pickle.dump(seq_motif_map, file, -1)
    assert os.path.isfile(output_path)
    if delete_intermediate:
        os.remove(conv_meme_file)
        shutil.rmtree(mast_meme_folder)
    elif storage_path:
        generic.warn_if_exist(storage_path, filetype='folder')
        if not os.path.isdir(storage_path):
            os.mkdir(storage_path)
        shutil.move(conv_meme_file, storage_path)
        shutil.move(mast_meme_folder, storage_path)

    return


#
#     assert isinstance(num_p, int)
#     assert num_p >= 1
#     assert process in ('meme', 'mast')
#     assert source in ('prosite', 'ioncom')
#
#     generic.quit_if_missing(paths.IONCOM_BINDING_SITES)
#     generic.warn_if_exist(paths.CONV_SEED_SEQS)
#     conv_interface.run(paths.IONCOM_BINDING_SITES, paths.CONV_SEED_SEQS)
#
#
#     parse_extracts(source, extract_path, pname_cid_path)
#     download_pdb(pname_cid_path, pdb_folder)
#     trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
#     create_seq(pname_cid_path, pdb_folder, seq_path)
#
#     filter_seq_file(seq_path)
#     # create_conv_seed_seqs()
#     find_motifs(process, pname_cid_path, ref_meme_txt, mast_meme_folder,
#                 seq_path, output, num_p)
#     if delete_intermediate:
#         os.remove(pname_cid_path)
#         os.remove(seq_path)
#         shutil.rmtree(mast_meme_folder)
#     elif storage_path:
#         generic.warn_if_exist(storage_path, filetype='folder')
#         if not os.path.isdir(storage_path):
#             os.mkdir(storage_path)
#         shutil.move(pname_cid_path, storage_path)
#         shutil.move(seq_path, storage_path)
#         shutil.move(mast_meme_folder, storage_path)
#
# def run_converge():
#     # Make seed_seq_file




def main():
    logs.set_logging_level()
    # Prosite, meme
    timestamp = datetime.datetime.now().isoformat()
    store_folder = os.path.join(paths.DEBUG, timestamp)
    os.mkdir(store_folder)
    output_path = os.path.join(paths.USER_OUTPUT, "prosite_meme_motif_pos.pkl")
    run_all(process='meme',
            source='prosite',
            num_p=7,
            extract_path=paths.PROSITE_EXTRACT,
            output=output_path,
            storage_path=store_folder)
    assert os.path.isfile(output_path)
    assert os.path.isdir(store_folder)

    # # Prosite, mast
    timestamp = datetime.datetime.now().isoformat()
    store_folder = os.path.join(paths.DEBUG, timestamp)
    os.mkdir(store_folder)
    output_path = os.path.join(paths.USER_OUTPUT, "prosite_mast_motif_pos.pkl")
    run_all(process='mast',
            source='prosite',
            extract_path=paths.PROSITE_EXTRACT,
            output=output_path,
            storage_path=store_folder)
    assert os.path.isfile(output_path)
    assert os.path.isdir(store_folder)

    # Ioncom, mast
    timestamp = datetime.datetime.now().isoformat()
    store_folder = os.path.join(paths.DEBUG, timestamp)
    os.mkdir(store_folder)
    output_path = os.path.join(paths.USER_OUTPUT, "ioncom_mast_motif_pos.pkl")
    run_all(process='mast',
            source='ioncom',
            extract_path=paths.IONCOM_EXTRACT,
            output=output_path,
            storage_path=store_folder)
    assert os.path.isfile(output_path)
    assert os.path.isdir(store_folder)


def run_all(process='meme', source='prosite', num_p=7, extract_path=None,
            pname_cid_path=paths.PNAME_CID, pdb_folder=paths.PDB_FOLDER,
            seq_path=paths.FULL_SEQS, ref_meme_txt=paths.REF_MEME_TXT,
            mast_meme_folder=paths.MEME_MAST_FOLDER, output=paths.MOTIF_POS,
            delete_intermediate=False, storage_path=None):
    assert isinstance(num_p, int)
    assert num_p >= 1
    assert process in ('meme', 'mast')
    assert source in ('prosite', 'ioncom')

    parse_extracts(source, extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, pdb_folder, seq_path)

    filter_seq_file(seq_path)
    # create_conv_seed_seqs()
    find_motifs(process, pname_cid_path, ref_meme_txt, mast_meme_folder,
                seq_path, output, num_p)
    if delete_intermediate:
        os.remove(pname_cid_path)
        os.remove(seq_path)
        shutil.rmtree(mast_meme_folder)
    elif storage_path:
        generic.warn_if_exist(storage_path, filetype='folder')
        if not os.path.isdir(storage_path):
            os.mkdir(storage_path)
        shutil.move(pname_cid_path, storage_path)
        shutil.move(seq_path, storage_path)
        shutil.move(mast_meme_folder, storage_path)


def parse_extracts(source='prosite', input_file_path=None,
                   pname_cid_path=paths.PNAME_CID):
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


def download_pdb(pname_cid_path=paths.PNAME_CID, pdb_folder=paths.PDB_FOLDER,
                 replace_existing=False):
    if replace_existing:
        logging.warning(f"Replacing .pdb files in PDB_folder <{pdb_folder}>. "
                        f"Downloading from scratch takes a while.")
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    if not os.path.isdir(pdb_folder):
        logging.warning(f"PDB_folder in <{pdb_folder}> not found. Downloading "
                        f"from scratch takes a while.")
        os.mkdir(pdb_folder)
    download_pdb_files.download(pname_cid_map, pdb_folder, False)


def trim_pnames_based_on_pdb(pname_cid_path=paths.PNAME_CID,
                             pdb_folder=paths.PDB_FOLDER):
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    download_pdb_files.trim_pname_cid(pname_cid_map, pdb_folder)
    _test_seq_cid_map(pname_cid_map)
    with open(pname_cid_path, 'wb') as file:
        pickle.dump(pname_cid_map, file, -1)
    generic.quit_if_missing(pname_cid_path)


def create_seq(pname_cid_path=paths.PNAME_CID,
               pdb_folder=paths.PDB_FOLDER, seq_path=paths.FULL_SEQS):
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    seqs = create_seq_file.extract_sequence(pdb_folder, pname_cid_map,
                                            AA3_to_AA1=generic.AA3_to_AA1)
    with open(seq_path, 'w') as file:
        file.writelines(seqs)
    create_seq_file.test_fasta_match_pdb(seq_path, pdb_folder,
                                         pname_cid_map, generic.AA3_to_AA1)


def filter_seq_file(seq_path=paths.FULL_SEQS):
    generic.quit_if_missing(seq_path)
    filter_seqs.delete_short_seqs(seq_path, threshold=30)


def find_motifs(process,
                pname_cid_path=paths.PNAME_CID,
                ref_meme_txt=paths.REF_MEME_TXT,
                mast_meme_folder=paths.MEME_MAST_FOLDER,
                seq_file=paths.FULL_SEQS,
                output=paths.MOTIF_POS,
                num_p=1):
    assert process in ('mast', 'meme')
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    motif_pos = motif_finder.find(pname_cid_map,
                                  process=process,
                                  num_p=num_p,
                                  ref_meme_txt=ref_meme_txt,
                                  mast_meme_folder=mast_meme_folder,
                                  seq_file=seq_file)
    generic.warn_if_exist(output)
    with open(output, 'wb') as file:
        pickle.dump(motif_pos, file, -1)


def create_conv_seed_seqs(binding_site_path=paths.IONCOM_BINDING_SITES,
                          seed_seq_path=paths.CONV_SEED_SEQS):
    generic.quit_if_missing(binding_site_path)
    generic.quit_if_missing(seed_seq_path)
    make_conv_seed_seqs.make(binding_site_path, seed_seq_path)


def _test_seq_cid_map(seq_cid_map):
    for pname, cid in seq_cid_map.items():
        assert isinstance(pname, str)
        assert isinstance(cid, str)
        assert len(pname) == 4
        assert len(cid) == 1
    return True

# main()
