import os
import numpy as np

from global_config import pdb_files_dir
from .parsers.loader import Loader
from .parsers.pdb_list_parser import parse_ptr_file

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
             ('LEU', 'L'),
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
    pdf_files = Loader(file_path_no_suffix)
    # hb2_files = Loader(file_path_no_suffix + ".hb2")

    ATOM = pdf_files.parse_with('ATOMParser')
    MODRES = pdf_files.parse_with('MODRESParser')
    HETATM = pdf_files.parse_with('HETATMParser')
    hb2  = pdf_files.parse_with('HbondParser')

    if not MODRES.empty:
        ATOM_res = ATOM.res.copy()
        ATOM.res = _MODRES_sub(ATOM_res, MODRES.res, MODRES.std_res_name)

    # ATOM.res, ATOM.sno = _inplace_AA3_substitution(
    # copy.deepcopy(ATOM.res), copy.deepcopy(ATOM.sno), AA3_TO_AA1)
    return ATOM, HETATM, hb2

def load_all(ptr_file):
    ptr_properties = parse_ptr_file(ptr_file)
    top_5 = list(ptr_properties.keys())[:10]
    ptr_properties2 = dict()
    for key in top_5:
        ptr_properties2[key] = ptr_properties[key]
    ptr_properties = ptr_properties2

    # print(ptr_properties)
    # import sys
    # sys.exit()
    return_data = dict()
    to_delete = []
    for p_name, properties in ptr_properties.items():
        print(p_name)
        try:
            filepath = os.path.join(pdb_files_dir, p_name)
            if 'sno_markers' not in properties:
                print(f"sno_markers not in properties, for file {p_name}.")
                raise AssertionError
            if not os.path.isfile(filepath):
                print(f"File {p_name} not found.")
                continue
            try:
                file_data = _load_data(filepath)
            except AssertionError as e:
                print(f"{e} | Flagged: {p_name}")
                continue
            ATOM, HETATM, hb2 = file_data
            if 'cid' in properties:
                ATOM = ATOM[ATOM.cid == properties['cid']]
                # HETATM = HETATM[HETATM.cid == properties['cid']]
            return_data[(p_name, tuple(properties['sno_markers']))] = (ATOM,
                                                                       HETATM, hb2)
        except:
            to_delete.append(p_name)
    print(f"<{to_delete}>")
    for i in to_delete:
        del ptr_properties[i]
        # print(properties['sno_markers'])
        # print(list(ATOM.sno))
        # print("end")
        #
        # for i in properties['sno_markers']:
        #     print(i)
        #     print(ATOM[ATOM.sno == i-1])
        #     print(list(ATOM[ATOM.sno == i-1].res)[0])
        #     print(list(ATOM[ATOM.sno == i].res)[0])
        #     print(list(ATOM[ATOM.sno == i+1].res)[0])
        #     print("")
        # print("\n")
    return return_data

def load_specific(ptr_file, specified_file):
    ptr_properties = parse_ptr_file(ptr_file)
    properties = ptr_properties[specified_file]
    assert 'sno_markers' in properties
    try:
        file_data = _load_data(os.path.join(pdb_files_dir, specified_file))
    except AssertionError as e:
        print(f"{e} | Flagged: {specified_file}")
        raise e
    ATOM, HETATM, hb2 = file_data
    if 'cid' in properties:
        ATOM = ATOM[ATOM.cid == properties['cid']]
        HETATM = HETATM[HETATM.cid == properties['cid']]
        hb2 = hb2[hb2.cid == properties['cid']]
    return_data = dict()
    return_data[(specified_file, properties['sno_markers'])] = \
        (ATOM, HETATM, hb2)
    return return_data

def load_all_after(ptr_file, specified_file):
    ptr_properties = parse_ptr_file(ptr_file)
    return_data = dict()
    start = False
    for p_name, properties in ptr_properties.items():
        if p_name == specified_file:
            start = True
        if start:
            file_data = _load_data(os.path.join(pdb_files_dir, p_name))
            ATOM, HETATM, hb2 = file_data
            if 'cid' in properties:
                ATOM = ATOM[ATOM.cid == properties['cid']]
                HETATM = HETATM[HETATM.cid == properties['cid']]
                hb2 = hb2[hb2.cid == properties['cid']]
            return_data[(p_name, properties['sno_markers'])] = \
                (ATOM, HETATM, hb2)
    return return_data
