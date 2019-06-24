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

Run setup_pname_cid_map_prosite/ioncom before
setup_motif_finder_prosite/ioncom, as motif_finder requires pname_cid_map.
"""
import os
import pickle
import shutil

import config
from utils import extract_parser, motif_finder

def check_if_input_files_exist():
    prosite_extract_path = os.path.join(
        config.ROOT, 'tests', 'data', 'input', 'prosite_extract.txt')
    ioncom_path = os.path.join(
        config.ROOT, 'tests', 'data', 'input', 'ioncom.txt')

    assert os.path.isfile(prosite_extract_path)
    assert os.path.isfile(ioncom_path)

def setup_pname_cid_map_prosite():
    prosite_extract_path = os.path.join(
        config.ROOT, 'tests', 'data', 'input', 'prosite_extract.txt')

    pname_cid_map_prosite = extract_parser.parse_prosite(
        prosite_extract_path, config.prosite_pdb_list)

    pname_cid_prosite_path = os.path.join(
        config.ROOT, 'tests', 'data', 'ref', 'pname_cid_prosite.pkl')

    with open(pname_cid_prosite_path, 'wb') as file:
        pickle.dump(pname_cid_map_prosite, file, -1)

def setup_pname_cid_map_ioncom():
    ioncom_path = os.path.join(
        config.ROOT, 'tests', 'data', 'input', 'ioncom.txt')

    pname_cid_map_ioncom = extract_parser.parse_ioncom(ioncom_path)

    pname_cid_ioncom_path = os.path.join(
        config.ROOT, 'tests', 'data', 'ref', 'pname_cid_ioncom.pkl')

    with open(pname_cid_ioncom_path, 'wb') as file:
        pickle.dump(pname_cid_map_ioncom, file, -1)

def setup_motif_finder_prosite():
    tmp_folder_path = os.path.join(config.ROOT, 'data', 'tmp')
    if os.path.isdir(tmp_folder_path):
        shutil.rmtree(tmp_folder_path)
    os.mkdir(tmp_folder_path)
    pdb_folder_path = os.path.join(config.ROOT, 'data', 'input',
                                   'pdb_files')

    pname_cid_prosite_path = os.path.join(
        config.ROOT, 'tests', 'data', 'ref', 'pname_cid_prosite.pkl')

    with open(pname_cid_prosite_path, 'rb') as file:
        pname_cid_prosite = pickle.load(file)

    ptr_props_prosite = motif_finder.find_motif_pos(
        pname_cid_prosite, pdb_folder_path, type='prosite',
        store_dir=tmp_folder_path, replace_existing=False,
        delete_intermediate_store=True, ref_meme_txt=None)

    ref_prosite_output_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                           'ptr_props_prosite.pkl')

    with open(ref_prosite_output_path, 'wb') as file:
        pickle.dump(ptr_props_prosite, file, -1)

def setup_motif_finder_ioncom():
    tmp_folder_path = os.path.join(config.ROOT, 'data', 'tmp')
    meme_txt_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                 'meme.txt')

    if os.path.isdir(tmp_folder_path):
        shutil.rmtree(tmp_folder_path)
    os.mkdir(tmp_folder_path)
    pdb_folder_path = os.path.join(config.ROOT, 'data', 'input',
                                   'pdb_files')

    pname_cid_ioncom_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                          'pname_cid_ioncom.pkl')

    with open(pname_cid_ioncom_path, 'rb') as file:
        pname_cid_ioncom = pickle.load(file)

    ptr_props_ioncom = motif_finder.find_motif_pos(
        pname_cid_ioncom,
        pdb_folder_path,
        type='ioncom',
        store_dir=tmp_folder_path,
        replace_existing=False,
        delete_intermediate_store=False,
        ref_meme_txt=meme_txt_path)

    ref_ioncom_output_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                           'ptr_props_ioncom.pkl')

    with open(ref_ioncom_output_path, 'wb') as file:
        pickle.dump(ptr_props_ioncom, file, -1)

if __name__ == '__main__':
    check_if_input_files_exist()
    setup_pname_cid_map_prosite()
    setup_pname_cid_map_ioncom()
    setup_motif_finder_prosite()
    setup_motif_finder_ioncom()
