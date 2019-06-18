import os
import logging
import numpy as np

import config
from descr.parsers import loader
from descr.parsers import pdb_list_parser

##############################################################################
# AA3_to_AA1
##############################################################################

# List of three and one letter amino acid codes
_aa_index = [('ALA', 'A'),
             ('CYS', 'C'),
             ('ASP', 'D'),
             ('GLU', 'E'),
             ('PHE', 'F'),
             ('GLY', 'G'),
             ('HIS', 'H'),
             ('HSE', 'H'),
             ('HSD', 'H'),
             ('ILE', 'I'),
             ('LYS', 'K'),
             ('LEU', 'L'),gl
             ('MET', 'M'),
             ('MSE', 'M'),
             ('ASN', 'N'),
             ('PRO', 'P'),
             ('GLN', 'Q'),
             ('ARG', 'R'),
             ('SER', 'S'),
             ('THR', 'T'),
             ('VAL', 'V'),
             ('TRP', 'W'),
             ('TYR', 'Y')]

AA3_TO_AA1 = dict(_aa_index)

def _inplace_AA3_substitution(res, sno, AA3_TO_AA1):
    """
    Not used, not deprecated.
    """
    assert len(res) == len(sno)

    for i in range(len(res)):
        try:
            res.iloc[i] = AA3_TO_AA1[res.iloc[i]]
        except KeyError:
            res.iloc[i] = "X"
            sno.iloc[i] = 0
    return res, sno

##############################################################################

def _MODRES_sub(main_res, MODRES_res, MODRES_std_res_name):
    for ATOM_i, ATOM_res in enumerate(main_res):
        for MODRES_i, res in enumerate(MODRES_res):
            if ATOM_res == res:
                main_res[ATOM_i] = MODRES_std_res_name[MODRES_i]
    return main_res

def _load_data(file_path_no_suffix):
    """
    For AA3_to_AA1, copy makes it much faster, because otherwise df will
    update itself with every change.
    """
    pdf_files = loader.Loader(file_path_no_suffix)
    # hb2_files = Loader(file_path_no_suffix + ".hb2")

    ATOM = pdf_files.parse_with('ATOMParser')
    MODRES = pdf_files.parse_with('MODRESParser')
    HETATM = pdf_files.parse_with('HETATMParser')
    hb2  = pdf_files.parse_with('HbondParser')

    if not MODRES.empty:
        ATOM_res = ATOM.res.copy()
        ATOM.res = _MODRES_sub(ATOM_res, MODRES.res, MODRES.std_res_name)

    return ATOM, HETATM, hb2

def load_pointer_file(ptr_file):
    # todo: fix parse_ptr_file so it either retrieve .pkl if exists, or rerun
    #  it from a main ptr script somewhere
    ptr_properties = pdb_list_parser.parse_ptr_file(ptr_file)
    # top_5 = list(ptr_properties.keys())[:10]
    # ptr_properties2 = dict()
    # for key in top_5:
    #     ptr_properties2[key] = ptr_properties[key]
    # ptr_properties = ptr_properties2

    return_data = dict()
    to_delete = []
    for i, (p_name, properties) in enumerate(ptr_properties.items()):
        print(f"{i}: {p_name}")
        try:
            filepath = os.path.join(config.pdb_files_dir, p_name+".pdb")
            if 'sno_markers' not in properties:
                logging.error(f"sno_markers not in properties, for "
                              f"file {p_name}.")
                raise AssertionError
            if not os.path.isfile(filepath):
                print(f"File {p_name} not found.")
                print(filepath)
                continue
            try:
                file_data = _load_data(filepath)
            except AssertionError as e:
                print(f"{e} | Flagged: {p_name}")
                continue
            ATOM, HETATM, hb2 = file_data
            # if 'cid' in properties:
            #     ATOM = ATOM[ATOM.cid == properties['cid']]
            return_data[(p_name, tuple(properties['sno_markers']),
                         properties['cid'])] = \
                (ATOM, HETATM, hb2)
        except Exception as e:
            print(e)
            print("kick")
            to_delete.append(p_name)
    print(f"<{to_delete}>")
    for i in to_delete:
        del ptr_properties[i]
    return return_data
