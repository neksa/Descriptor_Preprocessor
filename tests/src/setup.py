


# ----------------------------------------------
# Actual files setup

# We start off with the respective website extracts.
# For prosite, we need the pdb_list, obtained from copy-pasting from the
# website directly, and a html extract obtained by inspecting the page.
import os
import config
import pickle
import shutil

from utils import pdb_list, extract_parser

prosite_extract_path = os.path.join(config.ROOT, 'data', 'input',
                                    'prosite_extract.txt')

# Next, we convert both into a pname_cid_map, and save it in the input folder.
def setup_act_pname_cid_map_prosite():
    pname_cid_map_prosite = extract_parser.parse_prosite_extract(
        prosite_extract_path, config.prosite_pdb_list)

    pname_cid_prosite_path = os.path.join(config.ROOT, 'data', 'input',
                                        'pname_cid_prosite.pkl')

    with open(pname_cid_prosite_path, 'wb') as file:
        pickle.dump(pname_cid_map_prosite, file, -1)

# For datalon, the same procedure is performed, but we only need the website
# extract and not a pdb_list.
def setup_act_pname_cid_map_datalon():
    datalon_extract_path = os.path.join(config.ROOT, 'data', 'input',
                                        'allid_reso3.0_len50_nr40.txt')

    pname_cid_map_datalon = extract_parser.parse_datalon(datalon_extract_path)

    pname_cid_datalon_path = os.path.join(
        config.ROOT, 'tests', 'data', 'input',
                                        'pname_cid_datalon.pkl')

    with open(pname_cid_datalon_path, 'wb') as file:
        pickle.dump(pname_cid_map_datalon, file, -1)

# ----------------------------------------------
# Test Files Setup
#
# To obtain the reference output, we manually crop the website extracts to
# reduce runtime. We then perform the same extract procedure as before,
# saving the output file to the test ref folder.
#
# For prosite:

def setup_test_pname_cid_map_prosite():
    prosite_extract_shortened_path = os.path.join(
        config.ROOT, 'tests', 'data', 'input', 'prosite_extract_cropped.txt')

    pname_cid_map_prosite = extract_parser.parse_prosite_extract(prosite_extract_shortened_path,
                                                  config.prosite_pdb_list)

    pname_cid_prosite_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                        'pname_cid_prosite.pkl')

    with open(pname_cid_prosite_path, 'wb') as file:
        pickle.dump(pname_cid_map_prosite, file, -1)

# For datalon

def setup_test_pname_cid_map_datalon():
    datalon_shortened_path = os.path.join(
        config.ROOT, 'tests', 'data', 'input', 'datalon_shortened.txt')

    pname_cid_map_datalon = extract_parser.parse_datalon(datalon_shortened_path)

    pname_cid_datalon_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                        'pname_cid_datalon.pkl')

    with open(pname_cid_datalon_path, 'wb') as file:
        pickle.dump(pname_cid_map_datalon, file, -1)

# Then, we run get_pdb_list() to obtain a reference set of output, for both
# prosite and datalon. We use a meme.txt, obtained by running meme on the
# entire prosite extract previously (to obtain the motif pssm). The
# pdb_folder provided here is merely a tmp folder, to be deleted after. The
# output is saved as reference prosite output. Note that it takes a while for
# the pdb files to be downloaded. To short-circuit, use

# For prosite
def setup_test_ptr_props_prosite(reuse_pdb_folder=True):
    tmp_folder_path = os.path.join(config.ROOT, 'data', 'tmp')
    if os.path.isdir(tmp_folder_path):
        shutil.rmtree(tmp_folder_path)
    os.mkdir(tmp_folder_path)
    if reuse_pdb_folder:
        pdb_folder_path = os.path.join(config.ROOT, 'data', 'input', 'pdb_files')
    else:
        pdb_folder_path = os.path.join(config.ROOT, 'data', 'tmp', 'pdb')
        os.mkdir(pdb_folder_path)

    pname_cid_prosite_path = os.path.join(config.ROOT, 'tests', 'data',
                                          'ref',
                                          'pname_cid_prosite.pkl')

    with open(pname_cid_prosite_path, 'rb') as file:
        pname_cid_prosite = pickle.load(file)

    ptr_props_prosite = pdb_list.get_pdb_list(pname_cid_prosite,
                                              pdb_folder_path,
                                     type='prosite_extract',
                                     store_dir=tmp_folder_path,
                                     replace_existing=False,
                                     delete_intermediate_store=True,
                                     ref_meme_txt=None)

    ref_prosite_output_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                           'ptr_props_cropped_prosite.pkl')

    with open(ref_prosite_output_path, 'wb') as file:
        pickle.dump(ptr_props_prosite, file, -1)

    if not reuse_pdb_folder:
        shutil.rmtree(pdb_folder_path)

# For datalon
def setup_test_ptr_props_datalon(reuse_pdb_folder=True):
    tmp_folder_path = os.path.join(config.ROOT, 'data', 'tmp')
    meme_txt_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                 'meme.txt')

    if os.path.isdir(tmp_folder_path):
        shutil.rmtree(tmp_folder_path)
    os.mkdir(tmp_folder_path)
    if reuse_pdb_folder:
        pdb_folder_path = os.path.join(config.ROOT, 'data', 'input',
                                       'pdb_files')
    else:
        pdb_folder_path = os.path.join(config.ROOT, 'data', 'tmp', 'pdb')
        os.mkdir(pdb_folder_path)

    pname_cid_datalon_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                        'pname_cid_datalon.pkl')

    with open(pname_cid_datalon_path, 'rb') as file:
        pname_cid_datalon = pickle.load(file)

    ptr_props_datalon = pdb_list.get_pdb_list(pname_cid_datalon,
                                              pdb_folder_path,
                                              type='datalon',
                                              store_dir=tmp_folder_path,
                                              replace_existing=False,
                                              delete_intermediate_store=False,
                                              ref_meme_txt=meme_txt_path)

    ref_datalon_output_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                           'ptr_props_cropped_datalon.pkl')

    with open(ref_datalon_output_path, 'wb') as file:
        pickle.dump(ptr_props_datalon, file, -1)

    if not reuse_pdb_folder:
        shutil.rmtree(pdb_folder_path)

setup_test_pname_cid_map_prosite()
setup_test_pname_cid_map_datalon()
setup_test_ptr_props_prosite()
setup_test_ptr_props_datalon()

