import contextlib
import logging
import os
import pickle
import traceback
from urllib import request
import urllib.error

from pdb_component import pdb_utils, pdb_paths, loaders



def get_seq_for(pdb_code, cid=None):
    filedata = get_info_for(pdb_code)
    if filedata is None:
        return None
    ATOM, HETATM, hb = filedata
    del HETATM
    del hb
    if cid:
        ATOM_cid = ATOM[ATOM.cid == cid]
        seq = _extract_seq_from_df(ATOM_cid)
        return seq
    cid_seq_map = dict()
    for current_cid, ATOM_cid in ATOM.groupby("cid"):
        seq = _extract_seq_from_df(ATOM_cid)
        cid_seq_map[current_cid] = seq
    return cid_seq_map


def get_info_for(pdb_code):
    pdb_suffix = pdb_code.lower().strip() + ".pkl"
    if pdb_suffix not in pdb_paths.PDB_PARSED_SET:
        if pdb_code.lower().strip() + ".pdb" in pdb_paths.PDB_FILES_SET:
            get_success = loaders.load_pdb_info(pdb_code)
        else:
            get_success = download(pdb_code, silent=False)
            if get_success:
                pdb_paths.PDB_FILES_SET = set(os.listdir(pdb_paths.PDB_FILES))
                get_success = loaders.load_pdb_info(pdb_code)
        if not get_success:
            print(f"get_info_for(pdb_code) failed for {pdb_code}")
            return None
        pdb_paths.PDB_PARSED_SET = set(os.listdir(pdb_paths.PDB_PARSED))
    filepath = os.path.join(pdb_paths.PDB_PARSED, pdb_suffix)
    with open(filepath, 'rb') as file:
        output = pickle.load(file)
    return output


def preload_all():
    for filename in pdb_paths.PDB_FILES_SET:
        pdb_code = filename.split(".")[0]
        loaders.load_pdb_info(pdb_code)


def download(pdb_code, silent=False):
    pdb_code = pdb_code.lower().strip()

    url = pdb_utils.PDB_URL_TEMPLATE.format(pdb_code)
    print(url)
    output_path = os.path.join(pdb_paths.PDB_FILES, pdb_code+".pdb")
    try:
        with contextlib.closing(request.urlopen(url)) as contents:
            with open(output_path, 'w') as output_file:
                output_file.write(contents.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if not silent:
            logging.info(f"download() fails for file {output_path}. Probably "
                         f"invalid pdb_code.")
            logging.info(f"Traceback: <{traceback.format_exc()}>")
            logging.info(f"Error_msg: <{e}>\n")
        return False
    assert os.path.isfile(output_path)
    return True


def _extract_seq_from_df(df):
    # Assumption that df is screened for cid already, so res is unique
    seq = []
    snos = set(df.sno)
    max_sno = max(snos)
    current_sno = 1  # sno in df starts from 1
    while current_sno < max_sno:
        if current_sno in snos:
            AA3 = list(df[df.sno == current_sno].res)[0]
            try:
                AA1 = pdb_utils.AA3_to_AA1[AA3]
            except IndexError:
                AA1 = "X"
        else:
            AA1 = "X"
        current_sno += 1
        seq.append(AA1)
    seq = "".join(seq)
    return seq


