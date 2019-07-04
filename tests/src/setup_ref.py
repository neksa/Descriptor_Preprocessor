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

import config
from utils import extract_parser, motif_finder
from descr import loaders, descr_main

def setup_pname_cid_map_prosite(input_html_extract_path, pdb_list, output):
    pname_cid_map_prosite = extract_parser.parse_prosite(
        input_html_extract_path, pdb_list)
    with open(output, 'wb') as file:
        pickle.dump(pname_cid_map_prosite, file, -1)

def setup_pname_cid_map_ioncom(input_path, output):
    pname_cid_map_ioncom = extract_parser.parse_ioncom(input_path)
    with open(output, 'wb') as file:
        pickle.dump(pname_cid_map_ioncom, file, -1)

def setup_motif_finder_prosite(pname_cid_path, pdb_folder_path, output):
    tmp_folder_path = os.path.join(config.ROOT, 'data', 'tmp')
    if os.path.isdir(tmp_folder_path):
        shutil.rmtree(tmp_folder_path)
    os.mkdir(tmp_folder_path)

    with open(pname_cid_path, 'rb') as file:
        pname_cid_prosite = pickle.load(file)

    motif_pos_prosite = motif_finder.find_motif_pos(
        pname_cid_prosite, pdb_folder_path, process='meme',
        store_dir=tmp_folder_path, replace_existing=False,
        delete_intermediate_store=False, ref_meme_txt=None)

    with open(output, 'wb') as file:
        pickle.dump(motif_pos_prosite, file, -1)

def setup_motif_finder_ioncom(pname_cid_path, meme_txt_path, pdb_folder_path,
                              output):
    tmp_folder_path = os.path.join(config.ROOT, 'data', 'tmp')
    if os.path.isdir(tmp_folder_path):
        shutil.rmtree(tmp_folder_path)
    os.mkdir(tmp_folder_path)

    with open(pname_cid_path, 'rb') as file:
        pname_cid_ioncom = pickle.load(file)

    motif_pos_ioncom = motif_finder.find_motif_pos(
        pname_cid_ioncom,
        pdb_folder_path,
        process='mast',
        store_dir=tmp_folder_path,
        replace_existing=False,
        delete_intermediate_store=False,
        ref_meme_txt=meme_txt_path)

    with open(output, 'wb') as file:
        pickle.dump(motif_pos_ioncom, file, -1)

def setup_motif_finder_prosite_mast(pname_cid_path, meme_txt_path, pdb_folder_path,
                              output):
    tmp_folder_path = os.path.join(config.ROOT, 'data', 'tmp')
    if os.path.isdir(tmp_folder_path):
        shutil.rmtree(tmp_folder_path)
    os.mkdir(tmp_folder_path)

    with open(pname_cid_path, 'rb') as file:
        pname_cid_prosite = pickle.load(file)

    motif_pos_prosite = motif_finder.find_motif_pos(
        pname_cid_prosite, pdb_folder_path, process='mast',
        store_dir=tmp_folder_path, replace_existing=False,
        delete_intermediate_store=False, ref_meme_txt=meme_txt_path)

    with open(output, 'wb') as file:
        pickle.dump(motif_pos_prosite, file, -1)

def setup_descr(motif_pos_path, pdb_dir, output):
    assert os.path.isfile(motif_pos_path)
    with open(motif_pos_path, 'rb') as file:
        motif_pos = pickle.load(file)

    pdb_info_data_map = loaders.load_pdb_info(motif_pos, pdb_dir=pdb_dir)
    descrs = descr_main.calculate(pdb_info_data_map)

    with open(output, 'wb') as file:
        pickle.dump(descrs, file, -1)


def setup_all():
    input_ = os.path.join(config.ROOT, 'tests', 'data', 'input')
    ref_ = os.path.join(config.ROOT, 'tests', 'data', 'ref')

    prosite_extract_path = os.path.join(input_, 'prosite_extract.txt')
    pname_cid_prosite = os.path.join(ref_, 'pname_cid_prosite.pkl')

    ioncom_path = os.path.join(input_, 'ioncom.txt')
    pname_cid_ioncom = os.path.join(ref_, 'pname_cid_ioncom.pkl')
    # In main data folder, not in tests, can share, downloading takes a while.
    pdb_folder_path = os.path.join(config.ROOT, 'data', 'input', 'pdb_files')
    prosite_motif_pos = os.path.join(ref_, 'motif_pos_prosite.pkl')

    meme_txt_path = os.path.join(ref_, 'meme.txt')
    ioncom_motif_pos = os.path.join(ref_, 'motif_pos_ioncom.pkl')

    prosite_mast_motif_pos = os.path.join(ref_, 'motif_pos_meme_prosite.pkl')
    prosite_meme_descr = os.path.join(ref_, 'descr_prosite_meme.pkl')
    prosite_mast_descr = os.path.join(ref_, 'descr_prosite_mast.pkl')
    ioncom_mast_descr = os.path.join(ref_, 'descr_ioncom_mast.pkl')

    setup_pname_cid_map_prosite(prosite_extract_path, config.prosite_pdb_list,
                                pname_cid_prosite)
    setup_pname_cid_map_ioncom(ioncom_path, pname_cid_ioncom)
    setup_motif_finder_prosite(pname_cid_prosite, pdb_folder_path,
                               prosite_motif_pos)
    setup_motif_finder_ioncom(pname_cid_ioncom, meme_txt_path,
                              pdb_folder_path, ioncom_motif_pos)
    setup_motif_finder_prosite_mast(pname_cid_prosite, meme_txt_path,
                                    pdb_folder_path, prosite_mast_motif_pos)

    setup_descr(prosite_motif_pos, pdb_folder_path, prosite_meme_descr)
    setup_descr(prosite_mast_motif_pos, pdb_folder_path, prosite_mast_descr)
    setup_descr(ioncom_motif_pos, pdb_folder_path, ioncom_mast_descr)

if __name__ == '__main__':
    setup_all()
