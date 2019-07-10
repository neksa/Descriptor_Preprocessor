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
import shutil
import subprocess

import numpy as np
import matplotlib.pyplot as plt

from config import paths
from utils import seq_logo, generic

def find_motif_pos(seq_cid_map, pdb_folder, process='meme',
                   store_dir=paths.STORE, replace_existing=False,
                   delete_intermediate_store=False, ref_meme_txt=None):
    assert process in ('meme', 'mast')
    if process == 'mast':
        assert isinstance(ref_meme_txt, str)
        assert os.path.isfile(ref_meme_txt)
    directory = _create_internal_dir(store_dir, replace_existing)
    _test_seq_cid_map(seq_cid_map)
    generic.download_pdb_files(seq_cid_map, pdb_folder, replace_existing=False)
    _test_pdb_files_present(seq_cid_map, pdb_folder)
    seqs = _extract_sequence(pdb_folder, seq_cid_map)
    with open(directory['fasta_fpath'], 'w') as file:
        file.writelines(seqs)
    _delete_short_seqs(directory['fasta_fpath'], threshold=30)
    _test_fasta_match_pdb(directory['fasta_fpath'], pdb_folder, seq_cid_map)
    seq_motif_map = _build_seq_motif_map(process, directory, replace_existing,
                                         ref_meme_txt=ref_meme_txt)
    seq_motif_map = _adjust_motif_diagram(seq_motif_map)
    seq_motif_map = _delete_gapped_motifs(seq_motif_map,
                                          directory['fasta_fpath'])
    motif_pos = _assemble_motif_pos(seq_motif_map, seq_cid_map)
    if delete_intermediate_store:
        shutil.rmtree(store_dir)
    return motif_pos

def _create_internal_dir(store_dir, replace_existing):
    if not os.path.isdir(store_dir):
        os.mkdir(store_dir)
    elif replace_existing:
        shutil.rmtree(store_dir)
        os.mkdir(store_dir)
    directory = dict()
    fasta_fpath = os.path.join(store_dir, "prosite_seqs.fasta")
    meme_out = os.path.join(store_dir, "meme_out")
    mast_out = os.path.join(store_dir, "mast_out")

    directory['fasta_fpath'] = fasta_fpath
    directory['meme_out'] = meme_out
    directory['mast_out'] = mast_out
    return directory

#pylint: disable=invalid-name
def _extract_sequence(pdb_folder, pname_cid_map, AA3_to_AA1=generic.AA3_to_AA1):
    """
    Extracts the seq for each pname-cid from their .pdb file.
    """
    seqs = []
    for pname in pname_cid_map.keys():
        cid = pname_cid_map[pname]
        assert len(cid) == 1
        seq = [f">{pname}\n"]
        pdb_filepath = os.path.join(pdb_folder, pname+".pdb")
        with open(pdb_filepath, 'r') as pdb_file:
            current_atom_count = 1
            for line in pdb_file:
                if line.startswith("MODEL        2"):
                    break
                if not line.startswith("ATOM"):
                    continue
                if line[21] != cid:
                    continue
                if line[26] != " ": # altloc
                    continue
                sno = int(line[22:26].strip())
                if sno < current_atom_count:
                    continue
                elif sno > current_atom_count:
                    num_spacer_res = sno - current_atom_count
                    seq.append("X" * num_spacer_res)
                    current_atom_count = sno
                res = line[17:20].strip()
                AA_single_letter = AA3_to_AA1[res]
                seq.append(f"{AA_single_letter}")
                current_atom_count += 1
            seq.append("\n")
        seq = "".join(seq)
        seqs.append(seq)
    return seqs
# pylint: enable=invalid-name

def _delete_short_seqs(fasta_fname, threshold=30):
    """
    Delete sequences with len < threshold.
    # todo: change this s.t. it delete from data and not from fasta.
    """
    output_str = []
    tmp = ""
    with open(fasta_fname, 'r') as file:
        for line in file:
            if line.startswith(">"):
                tmp = line
                continue
            if len(line) < threshold:
                continue
            else:
                output_str.append(tmp + line)
                continue
    output_str = ''.join(output_str)
    with open(fasta_fname, 'w') as file:
        file.write(output_str)

def _run_meme(meme_output, fasta_filename, replace_existing=False):
    if replace_existing or not os.path.isdir(meme_output):
        command = f"meme -w 13 -protein -nmotifs 1 -mod anr " \
            f"-oc {meme_output} {fasta_filename}"
        subprocess.run(command, shell=True)
    return True

def _run_mast(mast_output, meme_txt, fasta_filename, replace_existing=False):
    if replace_existing or not os.path.isdir(mast_output):
        command = f"mast -oc {mast_output} -mt 0.0001 {meme_txt} " \
            f"{fasta_filename}"
        subprocess.run(command, shell=True)
    return True

def _get_motif_diagram(input_txt, source='meme'):
    """
    We obtain from meme/mast output (meme.txt/mast.txt) the motif diagram,
    showing the location of the matched motifs for each pname sequence. We
    then adjust it such that the locations are absolute, relative to the
    start of the sequence rather than from the end of the previous motif.
    Motifs with gaps in between are also deleted, mainly because the
    subsequent descriptor code assumes a continuous sequence for the len-30
    analysis segment.
    """
    seq_motif_map = dict()
    if source == 'meme':
        motif_diagram_sep = "_[1]_"
    elif source == "mast":
        motif_diagram_sep = "-[1]-"
    else:
        raise Exception
    with open(input_txt, 'r') as rfile:
        correct_segment = False
        in_area = False
        for line in rfile:
            if source == 'meme' and re.search("MEME-1 block diagrams", line):
                correct_segment = True
                continue
            elif source == 'mast' and line.startswith(
                    "SECTION II: MOTIF DIAGRAMS"):
                correct_segment = True
                continue
            if correct_segment and line.startswith("-------------   "):
                in_area = True
                continue
            if source == 'meme':
                if in_area and line.startswith("---------------------------"):
                    break
            else:
                if in_area and not line.strip():
                    break
            if in_area:
                pname, motif_diagram = line[:4], line[43:]
                if motif_diagram_sep not in motif_diagram:
                    continue
                raw_motif_positions = motif_diagram.split(motif_diagram_sep)
                motif_positions = []
                for pos in raw_motif_positions[:-1]:
                    if pos == "":
                        motif_positions.append(0)
                    else:
                        motif_positions.append(int(pos))
                seq_motif_map[pname] = motif_positions
    return seq_motif_map

def _adjust_motif_diagram(prev_map):
    seq_motif_map = dict()
    motif_len = 13
    for pname, relative_motif_pos in prev_map.items():
        abs_motif_pos = []
        curr_count = 1
        for i, pos in enumerate(relative_motif_pos):
            if i == 0:
                abs_pos = pos + 1
                curr_count += pos
                curr_count += motif_len
            else:
                abs_pos = curr_count + pos
                curr_count += pos
                curr_count += motif_len
            abs_motif_pos.append(abs_pos)
        seq_motif_map[pname] = abs_motif_pos
        assert len(relative_motif_pos) == len(abs_motif_pos)
    return seq_motif_map

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

import pickle
def _build_seq_motif_map(process, directory, replace_existing,
                         ref_meme_txt=None):
    if process == 'meme':
        _run_meme(directory['meme_out'], directory['fasta_fpath'],
                  replace_existing=replace_existing)
        _test_successful_meme(directory['meme_out'])
        meme_txt_path = os.path.join(directory['meme_out'], 'meme.txt')
        seq_motif_map = _get_motif_diagram(meme_txt_path, 'meme')
    elif process == 'mast':
        _run_mast(directory['mast_out'], ref_meme_txt,
                  directory['fasta_fpath'],
                  replace_existing=replace_existing)
        _test_successful_mast(directory['mast_out'])
        mast_txt_path = os.path.join(directory['mast_out'], 'mast.txt')
        seq_motif_map = _get_motif_diagram(mast_txt_path, 'mast')
    else:
        raise Exception
    return seq_motif_map

def _assemble_motif_pos(seq_motif_map, pname_cid_map):
    motif_positions = defaultdict(dict)
    for pname, each_pname_motif_pos in seq_motif_map.items():
        motif_positions[pname]['sno_markers'] = each_pname_motif_pos
        motif_positions[pname]['cid'] = pname_cid_map[pname]
    return motif_positions

def plot_labelled_motifs(seq_motif_map, fasta_filename):
    """
    This debugging function displays the sequence logo of the identified
    motifs.
    """
    motifs = _get_labelled_motifs(seq_motif_map, fasta_filename)
    plt.figure()
    converted_motifs = []
    extracted_motifs = np.array([list(i) for i in motifs])
    for AA_per_pos in extracted_motifs.T:
        AA, counts = np.unique(AA_per_pos, return_counts=True)
        sorted_i = np.argsort(counts)
        counts = counts[sorted_i]
        AA = AA[sorted_i]
        probs = counts / np.sum(counts)
        motif_this_pos = []
        for aa, prob in zip(AA, probs):
            motif_this_pos.append([aa, prob])
        converted_motifs.append(motif_this_pos)
    seq_logo.Logo(converted_motifs, -1, convert_AA3=False)

def _get_labelled_motifs(seq_motif_map, fasta_filename):
    labelled_motifs = []
    pname = None
    with open(fasta_filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                pname = line[1:5]
                continue
            if pname in seq_motif_map:
                motif_pos = seq_motif_map[pname]
            else:
                #   Might have been dropped because of gapped motif, etc
                pname = None
                continue
            for pos in motif_pos:
                if len(line) <= pos + 13:
                    logging.error(
                        f"Fasta seq is shorter than pos+13, for pos in "
                        f"motif_pos. Fasta_seq: <{line}>, "
                        f"motif_pos: <{motif_pos}>, illegal pos: <{pos}>.")
                    raise Exception
                labelled_motifs.append(line[pos - 1:pos + 12])
    return labelled_motifs

def _test_fasta_match_pdb(fasta_fname, pdb_folder, seq_cid_map,
                          AA3_to_AA1=generic.AA3_to_AA1):
    pname = None
    with open(fasta_fname, 'r') as f_file:
        for f_line in f_file:
            if f_line.startswith(">"):
                pname = f_line[1:5]
                continue
            if pname:
                pdb_path = os.path.join(pdb_folder, pname+'.pdb')
                cid = seq_cid_map[pname]
                with open(pdb_path, 'r') as p_file:
                    for p_line in p_file:
                        if not p_line.startswith("ATOM"):
                            continue
                        if p_line[21] != cid:
                            continue
                        if p_line[26] != " ":  # altloc
                            continue
                        sno = int(p_line[22:26].strip())
                        res = p_line[17:20]
                        if sno < 1:
                            continue
                        if AA3_to_AA1[res] != f_line[sno-1]:
                            logging.error(
                                f"Mismatch between seq in pdb file and in "
                                f"extracted fasta file. pdb_line: <{p_line}>, "
                                f"fasta_line: <{f_line}>")
                            raise Exception
    return True

def _test_seq_cid_map(seq_cid_map):
    for pname, cid in seq_cid_map.items():
        assert isinstance(pname, str)
        assert isinstance(cid, str)
        assert len(pname) == 4
        assert len(cid) == 1
    return True

def _test_pdb_files_present(seq_cid_map, pdb_folder, can_contain_extra=True):
    if not can_contain_extra:
        assert len(seq_cid_map) == len(os.listdir(pdb_folder))
    # check if all unique
    assert len(seq_cid_map) == len(set(seq_cid_map))

    fnames = []
    for fname in os.listdir(pdb_folder):
        split_fname = fname.split(".")
        assert len(split_fname) == 2
        pname = split_fname[0]
        fnames.append(pname)
    assert len(fnames) == len(set(fnames))
    fnames = set(fnames)
    for fname in seq_cid_map.keys():
        assert fname in fnames
    return True

def _test_successful_meme(meme_out):
    assert os.path.isdir(meme_out)
    _meme_txt_path = os.path.join(meme_out, 'meme.txt')
    assert os.path.isfile(_meme_txt_path)
    return True

def _test_successful_mast(mast_out):
    assert os.path.isdir(mast_out)
    _mast_txt_path = os.path.join(mast_out, 'mast.txt')
    assert os.path.isfile(_mast_txt_path)
    return True
