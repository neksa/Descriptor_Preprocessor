import logging
import os
import pickle
import shutil
from time import time

import config
from descr import descr_main, loaders
from utils import extract_parser, motif_finder, logs
import plots

def main():
    logs.set_logging_level()
    timecheck = time()

    # paths
    input_dir = os.path.join(config.ROOT, 'data', 'input')
    store_dir = os.path.join(config.ROOT, 'data', 'store')

    if os.path.isdir(store_dir):
        logging.warning("Store dir exists, deleting.")
        shutil.rmtree(store_dir)
    os.mkdir(store_dir)

    prosite_extract_path = os.path.join(input_dir, 'prosite',
                                        'prosite_extract.txt')
    ioncom_path = os.path.join(input_dir, 'ioncom','ioncom.txt')
    pdb_folder = os.path.join(input_dir, 'pdb_files')
    ref_meme_txt = os.path.join(input_dir, 'meme.txt')
    pname_cid_path = os.path.join(store_dir, 'pname_cid_map.pkl')
    motif_pos_path = os.path.join(store_dir, 'motif_pos.pkl')
    pdb_info_data_path = os.path.join(store_dir, 'pdb_info_data.pkl')

    # todo: set up local tmp folder, for calculate. Otherwise
    #  delete_intermediate_store now deletes my files. Or, somehow output it
    #  into output rather than keeping it in store.

    # todo: remember igor's point about re-generating the motifs from ioncom.

    # config
    source = 'ioncom'
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
        assert os.path.isfile(motif_pos_path)
        with open(motif_pos_path, 'rb') as file:
            motif_pos = pickle.load(file)
    else:
        if source == 'prosite':
            process = 'meme'
        elif source == 'ioncom':
            process = 'mast'
        else:
            raise Exception
        tmp_dir = os.path.join(config.ROOT, 'data', 'tmp')
        motif_pos = motif_finder.find_motif_pos(pname_cid_map,
                                                pdb_folder,
                                                process,
                                                replace_existing=True,
                                                delete_intermediate_store=False,
                                                store_dir=tmp_dir,
                                                ref_meme_txt=ref_meme_txt)
        if os.path.isfile(motif_pos_path):
            logging.warning(f"motif_pos in <{motif_pos_path}> exists. "
                            f"Replacing.")
        with open(motif_pos_path, 'wb') as file:
            pickle.dump(motif_pos, file, -1)

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
            pickle.dump(pdb_info_data_map, file, -1)

    descrs = descr_main.calculate(pdb_info_data_map)
    # for __, descr in descrs.groupby(['filename', 'cid', 'seq_marker']):
    #     calc_descr.write_descr(descr)

    logging.debug(time() - timecheck)

    # Switching back to pkl to avoid false float comparison failures.
    with open(os.path.join(store_dir, "descrs.pkl"), "wb") as pklfile:
        pickle.dump(descrs, pklfile, -1)

    plots.plot_signature_logo(descrs)



if __name__ == "__main__":
    main()