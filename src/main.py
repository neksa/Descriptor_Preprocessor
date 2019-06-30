import logging
import os
import pandas.testing
import pickle as pkl
from time import time

from descr import descr_main, loaders
import config
from utils import logs
import pickle
from utils import extract_parser, motif_finder

def main():
    logs.set_logging_level()
    timecheck = time()

    # paths
    input_ = os.path.join(config.ROOT, 'data', 'input')
    store_ = os.path.join(config.ROOT, 'data', 'store')

    prosite_extract_path = os.path.join(input_, 'prosite', 'prosite_extract.txt')
    ioncom_path = os.path.join(input_, 'ioncom','ioncom.txt')
    pdb_folder = os.path.join(input_, 'pdb_files')
    ref_meme_txt = os.path.join(input_, 'meme.txt')
    pname_cid_path = os.path.join(store_, 'pname_cid_map.pkl')
    motif_pos_path = os.path.join(store_, 'motif_pos.pkl')
    pdb_info_data_path = os.path.join(store_, 'pdb_info_data.pkl')

    # config
    source = 'prosite'
    to_load_pname_cid = False
    to_load_motif_pos = False
    to_load_pdb_info = False

    if to_load_pdb_info:
        to_load_pname_cid = True
        to_load_motif_pos = True
    if to_load_motif_pos:
        to_load_pname_cid = True
    assert source in ('prosite', 'ioncom')

    if to_load_pname_cid:
        assert os.path.isfile(pname_cid_path)
        with open(pname_cid_path, 'rb') as file:
            pname_cid_map = pickle.load(file)
    else:
        if source == 'prosite':
            pdb_list = config.prosite_pdb_list
            pname_cid_map = extract_parser.parse_prosite(prosite_extract_path, pdb_list)
        elif source == 'ioncom':
            pname_cid_map = extract_parser.parse_ioncom(ioncom_path)
        else:
            raise Exception
        if os.path.isfile(pname_cid_path):
            logging.warning(f"pname_cid in <{pname_cid_path}> exists. "
                            f"Replacing.")
        with open(pname_cid_path, 'wb') as file:
            pickle.dump(pname_cid_map, file, -1)

    if to_load_motif_pos:
        assert os.path.isfile(motif_pos_path)p
        with open(motif_pos_path, 'rb') as file:
            motif_pos = pickle.load(file)
    else:
        if source == 'prosite':
            process = 'meme'
        elif source == 'ioncom':
            process = 'mast'
        else:
            raise Exception
        motif_pos = motif_finder.find_motif_pos(pname_cid_map,
                                                pdb_folder,
                                                process,
                                                replace_existing=True,
                                                delete_intermediate_store=True,
                                                ref_meme_txt=ref_meme_txt)
        if os.path.isfile(motif_pos_path):
            logging.warning(f"motif_pos in <{motif_pos_path}> exists. "
                            f"Replacing.")
        with open(motif_pos_path, 'wb') as file:
            pickle.load(motif_pos, file, -1)

    if to_load_pdb_info:
        assert os.path.isfile(pdb_info_data_path)
        with open(pdb_info_data_path, 'rb') as file:
            pdb_info_data_map = pickle.load(file)
    else:
        pdb_info_data_map = loaders.load_pdb_info(motif_pos)
        if os.path.isfile(pdb_info_data_path):
            logging.warning(f"pdb_info_data in <{pdb_info_data_path}> exists. "
                            f"Replacing.")
        with open(pdb_info_data_path, 'wb') as file:
            pickle.load(pdb_info_data_map, file, -1)

    descrs = descr_main.calculate(pdb_info_data_map)
    # for __, descr in descrs.groupby(['filename', 'cid', 'seq_marker']):
    #     calc_descr.write_descr(descr)

    logging.debug(time() - timecheck)

    # Switching back to pkl to avoid false float comparison failures.
    with open(os.path.join(store_, "descrs.pkl"), "wb") as pklfile:
        pkl.dump(descrs, pklfile, -1)

    with open(os.path.join(store_, "reference.pkl"), "rb") as pklfile:
        reference = pkl.load(pklfile)

    pandas.testing.assert_frame_equal(descrs, reference)

if __name__ == "__main__":
    main()