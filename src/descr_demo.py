import logging
import os
from pandas.testing import assert_frame_equal
import pickle as pkl
from time import time

from descr.calc_descr import calc_descrs, write_descr
from global_config import store_dir, input_dir
from utils.logs import set_logging_level
from loaders import load_pointer_file
import pickle

if __name__ == "__main__":
    set_logging_level()

    timecheck = time()
    ptr_filename = os.path.join(input_dir, 'assembed_ptr_data.pkl')

    ptr_data = load_pointer_file(ptr_filename)

    with open("ptr_data.pkl", 'wb') as file:
        pickle.dump(ptr_data, file, -1)

    descrs = calc_descrs(ptr_data)
    # for __, descr in descrs.groupby(['filename', 'cid', 'seq_marker']):
    #     write_descr(descr)

    logging.debug(time() - timecheck)

    # Switching back to pkl to avoid false float comparison failures.
    with open(os.path.join(store_dir, "current2.pkl"), "wb") \
            as pklfile:
        pkl.dump(descrs, pklfile, protocol=-1)

    with open(os.path.join(store_dir, "reference.pkl"), "rb") \
            as pklfile:
        reference = pkl.load(pklfile)

    assert_frame_equal(descrs, reference)

    pass