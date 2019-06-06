from collections import namedtuple
import pickle

# namedtuple definition should be kept in module for pickle to work

import re

def parse_ptr_file(file_path):
    #         to_store[motif].append([p_name, cid, res_no])
    with open(file_path, 'rb') as file:
        ptr_file = pickle.load(file)
    return ptr_file

    for p_name, properties in ptr_file.items():

        cid = filename_cid_map[p_name]

        p_name = p_name.lower()
    #     displace, to match with the rest
        pdb_list.scop_class.append('scop_class')
        pdb_list.scop_fold.append('scop_fold')
        pdb_list.filenames.append(p_name,)
        pdb_list.cid.append(cid,)
        pdb_list.sno_markers.append(tuple(res_nos))
    #     todo: account for how some are shorter than 10.

    # with open(file_path) as file:
    #     lines = [line for line in file]
    # del lines[0]    # delete header
    # pdb_list = PDBList(*[[] for __ in range(4)])
    # for line in lines:
    #     scop_class, scop_fold, pdb_file_names, first_snos = line.split("\t", 3)
    #     first_snos = tuple([int(resi) for resi in first_snos.split(";")])
    #     pdb_list.scop_class.append(scop_class)
    #     pdb_list.scop_fold.append(scop_fold)
    #     pdb_list.filenames.append(pdb_file_names)
    #     pdb_list.sno_markers.append(first_snos)
    return pdb_list

# def parse_ptr_file(file_path):
#     #         to_store[motif].append([p_name, cid, res_no])
#     with open(file_path, 'rb') as file:
#         to_store = pickle.load(file)
#     pdb_list = PDBList(*[[] for __ in range(5)])
#
#     with open("../../filename_cid_map.pkl", 'rb') as file:
#         filename_cid_map = pickle.load(file)
#
#     for p_name, res_nos in to_store.items():
#         cid = filename_cid_map[p_name]
#
#         p_name = p_name.lower()
#     #     displace, to match with the rest
#         pdb_list.scop_class.append('scop_class')
#         pdb_list.scop_fold.append('scop_fold')
#         pdb_list.filenames.append(p_name,)
#         pdb_list.cid.append(cid,)
#         pdb_list.sno_markers.append(tuple(res_nos))
#     #     todo: account for how some are shorter than 10.
#
#     # with open(file_path) as file:
#     #     lines = [line for line in file]
#     # del lines[0]    # delete header
#     # pdb_list = PDBList(*[[] for __ in range(4)])
#     # for line in lines:
#     #     scop_class, scop_fold, pdb_file_names, first_snos = line.split("\t", 3)
#     #     first_snos = tuple([int(resi) for resi in first_snos.split(";")])
#     #     pdb_list.scop_class.append(scop_class)
#     #     pdb_list.scop_fold.append(scop_fold)
#     #     pdb_list.filenames.append(pdb_file_names)
#     #     pdb_list.sno_markers.append(first_snos)
#     return pdb_list