import logging
import os
import pandas.testing
import pickle as pkl
from time import time

from descr import calc_descr, loaders
import config
from utils import logs
import pickle

if __name__ == "__main__":
    logs.set_logging_level()

    timecheck = time()
    ptr_filename = os.path.join(config.input_dir, 'datalon.pkl')

    ptr_data = loaders.load_pointer_file(ptr_filename)

    with open("ptr_data2.pkl", 'wb') as file:
        pickle.dump(ptr_data, file, -1)

    descrs = calc_descr.calc_descrs(ptr_data)
    # for __, descr in descrs.groupby(['filename', 'cid', 'seq_marker']):
    #     calc_descr.write_descr(descr)

    logging.debug(time() - timecheck)

    # Switching back to pkl to avoid false float comparison failures.
    with open(os.path.join(config.store_dir, "current2.pkl"), "wb") \
            as pklfile:
        pkl.dump(descrs, pklfile, protocol=-1)

    with open(os.path.join(config.store_dir, "reference.pkl"), "rb") \
            as pklfile:
        reference = pkl.load(pklfile)

    pandas.testing.assert_frame_equal(descrs, reference)

    pass