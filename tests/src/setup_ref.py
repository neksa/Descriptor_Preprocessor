"""
Produces the reference output files for test_motif_finder.py.

We have as input some form of information detailing the pname => cid mapping
and we wish to obtain a pname+cid => motif_pos at the end. First, we parse
the input information to obtain pname_cid_map. Together with a pre-downloaded
list of pdb files in pdb_folder, and a pre-computed meme.txt obtained from
running meme on the full ef-hand dataset or obtained elsewhere, this is fed
into parse_ioncom() to obtain the desired output. This is then saved
separately for each source.

For prosite, the list of pdb codes representing the sequences are in
https://prosite.expasy.org/PS00018, copying the >>more list in PDB tab,
in the table. This forms prosite_pdb_list in config. We then take the html by
view_page_source in
https://prosite.expasy.org/cgi-bin/pdb/pdb_structure_list.cgi?src=PS00018,
and copy it into prosite_extract.txt.

For ioncom, we download the dataset from:
https://zhanglab.ccmb.med.umich.edu/IonCom/. This already provides us with
the pname+cid directly, so there's no need for a separate pdb_list.

Run setup_pname_cid_map_prosite/ioncom() before
setup_motif_finder_prosite/ioncom(), as motif_finder requires pname_cid_map.

Run plot_ref.py after to visualise the seq logo of the derived descr.
"""
import os
import pickle
import shutil

from config import paths
from descr import loaders, descr_main

from tests.src import paths_test

from preprocessing import preprocess

# Preprocessing
def setup_preprocessing():
    setup_parse_extracts()
    setup_create_seq()
    setup_find_motif()
    setup_run_all()

def setup_parse_extracts():
    prosite_input = paths_test.PROSITE_EXTRACT
    prosite_ref = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
    preprocess.parse_extracts('prosite', prosite_input, prosite_ref)
    assert os.path.isfile(prosite_ref)

    ioncom_input = paths_test.IONCOM_EXTRACT
    ioncom_ref = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
    preprocess.parse_extracts('ioncom', ioncom_input, ioncom_ref)
    assert os.path.isfile(ioncom_ref)


def setup_create_seq():
    input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
    input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
    seq_1 = paths_test.REF_CREATE_SEQ_1
    seq_2 = paths_test.REF_CREATE_SEQ_2
    preprocess.create_seq(pname_cid_path=input_1,
                          pdb_folder=paths_test.PDB_FOLDER,
                          seq_path=seq_1)
    preprocess.create_seq(pname_cid_path=input_2,
                          pdb_folder=paths_test.PDB_FOLDER,
                          seq_path=seq_2)
    assert os.path.isfile(seq_1)
    assert os.path.isfile(seq_2)


def setup_find_motif():
    input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
    input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
    seq_1 = paths_test.REF_CREATE_SEQ_1
    seq_2 = paths_test.REF_CREATE_SEQ_2
    output_1 = paths_test.REF_FIND_MOTIF_1
    output_2 = paths_test.REF_FIND_MOTIF_2
    preprocess.find_motifs('meme',
                           pname_cid_path=input_1,
                           ref_meme_txt=None,
                           seq_file=seq_1,
                           output=output_1,
                           num_p=7)
    preprocess.find_motifs('mast',
                           pname_cid_path=input_2,
                           ref_meme_txt=paths_test.REF_MEME_TXT,
                           seq_file=seq_2,
                           output=output_2,
                           num_p=7)
    assert os.path.isfile(output_1)
    assert os.path.isfile(output_2)


def setup_run_all():
    input_1 = paths_test.PROSITE_EXTRACT
    input_2 = paths_test.IONCOM_EXTRACT
    output_1 = paths_test.REF_RUN_ALL_1
    output_2 = paths_test.REF_RUN_ALL_2
    output_3 = paths_test.REF_RUN_ALL_3
    ref_meme_txt = paths_test.REF_MEME_TXT
    preprocess.run_all(process='meme',
                       source='prosite',
                       num_p=7,
                       extract_path=input_1,
                       output=output_1)
    preprocess.run_all(process='mast',
                       source='prosite',
                       extract_path=input_1,
                       ref_meme_txt=ref_meme_txt,
                       output=output_2)
    preprocess.run_all(process='mast',
                       source='ioncom',
                       extract_path=input_2,
                       ref_meme_txt=ref_meme_txt,
                       output=output_3)
    assert os.path.isfile(output_1)
    assert os.path.isfile(output_2)
    assert os.path.isfile(output_3)

if __name__ == "__main__":
    setup_preprocessing()

# setup_parse_extracts()
# setup_create_seq()
# def setup_run_all():
#     pdb_folder = paths_test.PDB_FOLDER
#     seq_path = paths_test.FULL_SEQS
#     ref_meme_txt = None
#
#
#     preprocess.run_all(process='meme', num_p=7, extract_path=None,
#             pname_cid_path=paths.PNAME_CID, pdb_folder=paths.PDB_FOLDER,
#             seq_path=paths.FULL_SEQS, ref_meme_txt=paths.REF_MEME_TXT,
#             mast_meme_folder=paths.MEME_MAST_FOLDER)
#
# def setup_pname_cid_map_ioncom(input_path, output):
#     pname_cid_map_ioncom = extract_parser.parse_ioncom(input_path)
#     with open(output, 'wb') as file:
#         pickle.dump(pname_cid_map_ioncom, file, -1)
#
# def setup_motif_finder_prosite(pname_cid_path, pdb_folder_path, output):
#     tmp_folder_path = os.path.join(paths.ROOT, 'data', 'tmp')
#     if os.path.isdir(tmp_folder_path):
#         shutil.rmtree(tmp_folder_path)
#     os.mkdir(tmp_folder_path)
#
#     with open(pname_cid_path, 'rb') as file:
#         pname_cid_prosite = pickle.load(file)
#
#     motif_pos_prosite = motif_finder.find_motif_pos(
#         pname_cid_prosite, pdb_folder_path, process='meme',
#         store_dir=tmp_folder_path, replace_existing=False,
#         delete_intermediate_store=False, ref_meme_txt=None)
#
#     with open(output, 'wb') as file:
#         pickle.dump(motif_pos_prosite, file, -1)
#
# def setup_motif_finder_ioncom(pname_cid_path, meme_txt_path, pdb_folder_path,
#                               output):
#     tmp_folder_path = os.path.join(paths.ROOT, 'data', 'tmp')
#     if os.path.isdir(tmp_folder_path):
#         shutil.rmtree(tmp_folder_path)
#     os.mkdir(tmp_folder_path)
#
#     with open(pname_cid_path, 'rb') as file:
#         pname_cid_ioncom = pickle.load(file)
#
#     motif_pos_ioncom = motif_finder.find_motif_pos(
#         pname_cid_ioncom,
#         pdb_folder_path,
#         process='mast',
#         store_dir=tmp_folder_path,
#         replace_existing=False,
#         delete_intermediate_store=False,
#         ref_meme_txt=meme_txt_path)
#
#     with open(output, 'wb') as file:
#         pickle.dump(motif_pos_ioncom, file, -1)
#
# def setup_motif_finder_prosite_mast(pname_cid_path, meme_txt_path, pdb_folder_path,
#                               output):
#     tmp_folder_path = os.path.join(paths.ROOT, 'data', 'tmp')
#     if os.path.isdir(tmp_folder_path):
#         shutil.rmtree(tmp_folder_path)
#     os.mkdir(tmp_folder_path)
#
#     with open(pname_cid_path, 'rb') as file:
#         pname_cid_prosite = pickle.load(file)
#
#     motif_pos_prosite = motif_finder.find_motif_pos(
#         pname_cid_prosite, pdb_folder_path, process='mast',
#         store_dir=tmp_folder_path, replace_existing=False,
#         delete_intermediate_store=False, ref_meme_txt=meme_txt_path)
#
#     with open(output, 'wb') as file:
#         pickle.dump(motif_pos_prosite, file, -1)
#
# def setup_descr(motif_pos_path, pdb_dir, output):
#     assert os.path.isfile(motif_pos_path)
#     with open(motif_pos_path, 'rb') as file:
#         motif_pos = pickle.load(file)
#
#     pdb_info_data_map = loaders.load_pdb_info(motif_pos, pdb_dir=pdb_dir)
#     descrs = descr_main.calculate(pdb_info_data_map)
#
#     with open(output, 'wb') as file:
#         pickle.dump(descrs, file, -1)
#
#
# def setup_all():
#     input_ = os.path.join(paths.ROOT, 'tests', 'data', 'input')
#     ref_ = os.path.join(paths.ROOT, 'tests', 'data', 'ref')
#
#     prosite_extract_path = os.path.join(input_, 'prosite_extract.txt')
#     pname_cid_prosite = os.path.join(ref_, 'pname_cid_prosite.pkl')
#
#     ioncom_path = os.path.join(input_, 'ioncom.txt')
#     pname_cid_ioncom = os.path.join(ref_, 'pname_cid_ioncom.pkl')
#     # In main data folder, not in tests, can share, downloading takes a while.
#     pdb_folder_path = os.path.join(paths.ROOT, 'data', 'input', 'pdb_files')
#     prosite_motif_pos = os.path.join(ref_, 'motif_pos_prosite.pkl')
#
#     meme_txt_path = os.path.join(ref_, 'meme.txt')
#     ioncom_motif_pos = os.path.join(ref_, 'motif_pos_ioncom.pkl')
#
#     prosite_mast_motif_pos = os.path.join(ref_, 'motif_pos_meme_prosite.pkl')
#     prosite_meme_descr = os.path.join(ref_, 'descr_prosite_meme.pkl')
#     prosite_mast_descr = os.path.join(ref_, 'descr_prosite_mast.pkl')
#     ioncom_mast_descr = os.path.join(ref_, 'descr_ioncom_mast.pkl')
#
#     setup_pname_cid_map_prosite(prosite_extract_path, generic.prosite_pdb_list,
#                                 pname_cid_prosite)
#     setup_pname_cid_map_ioncom(ioncom_path, pname_cid_ioncom)
#     setup_motif_finder_prosite(pname_cid_prosite, pdb_folder_path,
#                                prosite_motif_pos)
#     setup_motif_finder_ioncom(pname_cid_ioncom, meme_txt_path,
#                               pdb_folder_path, ioncom_motif_pos)
#     setup_motif_finder_prosite_mast(pname_cid_prosite, meme_txt_path,
#                                     pdb_folder_path, prosite_mast_motif_pos)
#
#     setup_descr(prosite_motif_pos, pdb_folder_path, prosite_meme_descr)
#     setup_descr(prosite_mast_motif_pos, pdb_folder_path, prosite_mast_descr)
#     setup_descr(ioncom_motif_pos, pdb_folder_path, ioncom_mast_descr)
#
# if __name__ == '__main__':
#     setup_all()
