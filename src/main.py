import logging
import os
import pickle
import shutil
from time import time

from config import paths
from descr import descr_main
from utils import logs, plots, generic

def main():
    logs.set_logging_level()

    # paths
    store_dir = os.path.join(paths.ROOT, 'data', 'store')

    if os.path.isdir(store_dir):
        logging.warning("Store dir exists, deleting.")
        shutil.rmtree(store_dir)
    os.mkdir(store_dir)

    pdb_info_data_path = os.path.join(store_dir, 'pdb_info_data.pkl')

    generic.quit_if_missing(paths.OUTPUT_DESCRS)
    with open(pdb_info_data_path, 'rb') as file:
        pdb_info_data_map = pickle.load(file)

    timecheck = time()
    descrs = descr_main.calculate(pdb_info_data_map)
    logging.debug(f"Time taken: {time() - timecheck}")
    # for __, descr in descrs.groupby(['filename', 'cid', 'seq_marker']):
    #     calc_descr.write_descr(descr)


    generic.warn_if_exist(paths.OUTPUT_DESCRS)
    # Switching back to pkl to avoid false float comparison failures.
    with open(paths.OUTPUT_DESCRS, "wb") as file:
        pickle.dump(descrs, file, -1)

    plots.plot_signature_logo(descrs)



if __name__ == "__main__":
    main()