"""
Find the starting position of a motif, given either a set of matched
sequences or a pre-computed motif file.

If we have to derive the motif, meme/converge is run on the set of seq+cid
for which we know a match exists. Otherwise, if given a motif file (in
meme.txt), we run mast on the set of seq+cid. Either way, we get at the end
a map between seq+cid, and the starting position(s) of the motif.

If the pdb files do not exist, they're downloaded from rscb. This takes a
while, so the pdb files are shared and persisted as much as possible
throughout the src.

Additional:
1. sno positions should only be displayed for motifs that have a continuous
representation of offset[0] => offset[1]
2. sno positions should only be displayed for the relevant cid => corollary
that we either only do one cid per pname, or that it must be separately
listed and labelled (e.g. pname_A, pname_b,...)
"""
from collections import defaultdict
import logging
import os
import re
import subprocess

from utils import generic
from meme_suite import meme_interface

def find(pname_cid_map, motif_len, process,
         num_p, ref_meme_txt, meme_folder, seq_file):
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    assert process in ('meme', 'mast')
    if process == 'mast':
        assert isinstance(ref_meme_txt, str)
        assert os.path.isfile(ref_meme_txt)
    if process == 'meme':
        meme_interface.run_meme(seq_file, motif_len, meme_folder, num_p)
        meme_txt_path = os.path.join(meme_folder, 'meme.txt')
        seq_motif_map = meme_interface.extract_motifs_meme(meme_txt_path,
                                                           motif_len)
    else:
        meme_interface.run_mast(ref_meme_txt, seq_file, meme_folder)
        mast_txt_path = os.path.join(meme_folder, 'mast.txt')
        seq_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
                                                           motif_len)
    seq_motif_map = _delete_gapped_motifs(seq_motif_map, seq_file)
    motif_pos = _assemble_motif_pos(seq_motif_map, pname_cid_map)
    return motif_pos

def _delete_gapped_motifs(prev_map, fasta_fname):
    seq_motif_map = dict()
    pname = None
    with open(fasta_fname, 'r') as file:
        for line in file:
            if line.startswith(">"):
                pname = line[1:].strip()
                continue
            if pname and pname in prev_map:
                screened_motif_pos = []
                motif_pos = prev_map[pname]
                for pos in motif_pos:
                    try:
                        motif = line[pos-8:pos+22]
                    except IndexError:
                        logging.info(
                            f"{pos} in {pname} has IndexError when obtaining "
                            f"motif from -8 to 22, skipping.\n")
                        continue
                    else:
                        if "X" not in motif:    # It's a continuous seq
                            screened_motif_pos.append(pos)
                if screened_motif_pos:
                    seq_motif_map[pname] = screened_motif_pos
    return seq_motif_map

def _build_seq_motif_map(process, tmp_output_folder, seq_file, motif_len,
                         num_p=1,
                         ref_meme_txt=None):
    if process == 'meme':
        meme_interface.run_meme(seq_file, tmp_output_folder, num_p)
        meme_txt_path = os.path.join(tmp_output_folder, 'meme.txt')
        seq_motif_map = meme_interface.extract_motifs_meme(meme_txt_path,
                                                           motif_len)
    elif process == 'mast':
        meme_interface.run_mast(ref_meme_txt, seq_file, tmp_output_folder)
        mast_txt_path = os.path.join(tmp_output_folder, 'mast.txt')
        seq_motif_map = meme_interface.extract_motifs_mast(mast_txt_path,
                                                           motif_len)
    else:
        raise Exception
    return seq_motif_map

def _assemble_motif_pos(seq_motif_map, pname_cid_map):
    motif_positions = defaultdict(dict)
    for pname, each_pname_motif_pos in seq_motif_map.items():
        motif_positions[pname]['sno_markers'] = each_pname_motif_pos
        motif_positions[pname]['cid'] = pname_cid_map[pname]
    return motif_positions
