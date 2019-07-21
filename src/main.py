import logging
import os
import pickle
import shutil
from time import time

from config import paths
from descr import descr_main, loaders
from utils import logs, generic, plots

def main():
    logs.set_logging_level()
    timecheck = time()

    # paths
    store_dir = os.path.join(paths.ROOT, 'data', 'store')

    if os.path.isdir(store_dir):
        logging.warning("Store dir exists, deleting.")
        shutil.rmtree(store_dir)
    os.mkdir(store_dir)

    pdb_info_data_path = os.path.join(store_dir, 'pdb_info_data.pkl')

    with open(pdb_info_data_path, 'rb') as file:
        pdb_info_data_map = pickle.load(file)

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