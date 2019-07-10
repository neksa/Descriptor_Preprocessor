
import pickle
import os
import logging


from config import paths
from utils import generic
from preprocessing import extract_parser

def parse_extracts(source='prosite', filename=None):
    assert source in ('prosite', 'ioncom')

    if filename == None:
        if source == 'prosite':
            extract_path = paths.PROSITE_EXTRACT_PATH
        elif source == 'ioncom':
            extract_path = paths.IONCOM_EXTRACT_PATH
    else:
        extract_path = os.path.join(paths.USER_INPUT, filename)
    # noinspection PyUnboundLocalVariable
    if not os.path.isfile(extract_path):
        logging.error(f"Input extract file not found in <{extract_path}>. "
                      f"Exiting.")
        raise Exception

    if source == 'prosite':
        pname_cid_map = extract_parser.parse_prosite(extract_path,
                                                     generic.prosite_pdb_list)
    elif source == 'ioncom':
        pname_cid_map = extract_parser.parse_ioncom(extract_path)

    if os.path.isfile(paths.PNAME_CID):
        logging.warning(f"pname_cid in <{paths.PNAME_CID}> exists. "
                        f"Replacing.")
    with open(paths.PNAME_CID, 'wb') as file:
        # noinspection PyUnboundLocalVariable
        pickle.dump(pname_cid_map, file, -1)

    return