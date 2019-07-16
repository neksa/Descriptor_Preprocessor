"""
Parsers to extract pname_cid_map from various sources

For prosite, we take from website, both copy-pasting of pdb_list, and html
from page source.

For ioncom, we download and parse the sequence-cid list directly.
"""

import logging
import re

def parse_prosite(ptr_file, pdb_list):
    with open(ptr_file, 'r') as file:
        cids = []
        for line in file:
            if re.search("MOL_ID: 1", line):
                descr_for_id_1 = line.split("MOL_ID: 1")[1]
                value = re.search("CHAIN: ([A-Z1-9])", descr_for_id_1)
                if value is None:
                    logging.warning(f"Valid cid not found in prosite_extract "
                                    f"for line <{line}>.")
                cids.append(value[1])
    pname_cid_map = dict()
    for pname, cid in zip(pdb_list, cids):
        pname = pname.lower()
        pname_cid_map[pname] = cid
    return pname_cid_map

def parse_ioncom(ptr_file):
    # Example: 1abqA\n => pname, cid = 1abq, A
    pname_cid_map = dict()
    with open(ptr_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                assert len(line) == 5
                pname, cid = line[:4], line[4]
                pname_cid_map[pname] = cid
    return pname_cid_map
