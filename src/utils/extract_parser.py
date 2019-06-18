import re
import logging

def parse_prosite_extract(ptr_file, pdb_list):
    with open(ptr_file, 'r') as file:
        cids = []
        for line in file:
            if re.search("MOL_ID: 1", line):
                start, remain = line.split("MOL_ID: 1")
                value = re.search("CHAIN: ([A-Z1-9])", remain)
                if value is None:
                    logging.warning(f"Valid cid not found in prosite_extract "
                                    f"for line <{line}>.")
                cids.append(value.group(1))
    pname_cid_map = dict()
    for pname, cid in zip(pdb_list, cids):
        pname = pname.lower()
        pname_cid_map[pname] = cid
    return pname_cid_map

# ----------------------------------
# We take as input the datalon pdb list, and obtain from it the mapping between
# pname and cid to be used.
#
# Example: 1abqA\n => pname, cid = line[:4], line[4]

def parse_datalon(ptr_file):
    pname_cid_map = dict()
    with open(ptr_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                assert len(line) == 5
                pname, cid = line[:4], line[4]
                pname_cid_map[pname] = cid
    return pname_cid_map
