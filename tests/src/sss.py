import os
import config
import pickle
import sys

# datalon_input_path = os.path.join(config.ROOT, 'tests', 'data',
#                                           'input', 'datalon_test.txt')
#
# from utils.extract_parser import parse_prosite_extract, parse_datalon
#
# pname_cid_map_datalon = parse_datalon(datalon_input_path)
pname_cid_path = os.path.join(config.ROOT, 'tests', 'data',
                              'input', 'pname_cid_datalon.pkl')
#
# with open(pname_cid_path, 'wb') as file:
#     pickle.dump(pname_cid_map_datalon, file, -1)

with open(pname_cid_path, 'rb') as file:
    pname_cid_map = pickle.load(file)

ptr_props_path = os.path.join(config.ROOT, 'tests', 'data',
                              'ref', 'ptr_props.pkl')
from utils import pdb_list

pdb_folder = "tmp"

print(pname_cid_map)

ptr_props = pdb_list.get_pdb_list(
            pname_cid_map, pdb_folder, type='prosite_extract',
            store_dir=config.store_dir, replace_existing=True)
print(ptr_props)

with open(ptr_props, 'wb') as file:
    pickle.dump(ptr_props, file, -1)
    # pname_cid_map = pickle.load(file)

sys.exit()


from utils import pdb_list
# from config import pdb_list

prosite_extract_path = os.path.join(config.ROOT, 'data', 'input',
                                    'prosite_extract.txt')

datalon_input_path = os.path.join(config.ROOT, 'data', 'input',
                                    'allid_reso3.0_len50_nr40.txt')

# pname_cid_map_prosite = parse_prosite_extract(prosite_extract_path, pdb_list)

pname_cid_map_datalon = parse_datalon(datalon_input_path)

pname_cid_path = os.path.join(config.ROOT, 'tests', 'data',
                              'input', 'pname_cid_map_datalon.pkl')
# with open(pname_cid_path, 'rb') as file:
#     pname_cid_map = pickle.load(file)
pdb_folder = os.path.join(config.store_dir, 'pdb_folder')
print(pname_cid_map_datalon)
# ptr_props = pdb_list.get_pdb_list(
#             pname_cid_map, pdb_folder, type='prosite_extract',
#             store_dir=config.store_dir, replace_existing=True)
# print(ptr_props)

# data = parse_datalon(datalon_input_path)
# with open("ptr_props.pkl", 'wb') as file:
#     pickle.dump(ptr_props, file, -1)
with open(pname_cid_path, 'wb') as file:
    pickle.dump(pname_cid_map_datalon, file, -1)
