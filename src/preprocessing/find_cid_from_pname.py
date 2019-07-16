"""
Given a pdb_code:seq map and .pdb files, find corresponding cid for each code.

Does so by matching sequence to that in .pdb file.
"""
import os

import leven
import numpy as np

from config import paths

def find(pname_seq_map, pdb_folder=paths.PDB_FOLDER):
    # pname is pdb_code
    # all pnames should be present in pdb_folder by this stage, do a check
# before calling this func.
    pname_cid_map = dict()
    for pname, seq in pname_seq_map.items():
        pdb_filepath = os.path.join(pdb_folder, pname+".pdb")
        cid_seq_map = _extract_seq_from_pdb(pdb_filepath)
        cids = list(cid_seq_map.keys())
        assert len(cids) >= 1
        match_scores = []
        for compared_seq in cid_seq_map.values():
            score = _get_seq_match_score(seq, compared_seq)
            match_scores.append(score)
        match_i = np.argmin(match_scores)
        cid = cids[match_i]
        pname_cid_map[pname] = cid
    return pname_cid_map

def _get_seq_match_score(seq_1, seq_2):
    """
    Can use any reasonable matching criteria. Using leven for now, rather
    overkill for what we need but it's k for now.
    """
    return leven.levenshtein(seq_1, seq_2)



def _extract_seq_from_pdb(pdb_filepath):
    return None