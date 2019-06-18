"""
get_ptr_list(*args) should produce a dict, with keys as filenames of pdb files,
and values as a list of sno positions from which the pattern starts, without
offset.

args:
pdb_folder = location of pdb files
extract_file = location of ptr file to be parsed
type = ("prosite_extract",)

Additional:
1. sno positions should only be displayed for motifs that have a continuous
representation of offset[0] => offset[1]
2. sno positions should only be displayed for the relevant cid => corollary
that we either only do one cid per pname, or that it must be separately
listed and labelled (e.g. pname_A, pname_b,...)

"""
# todo: change Logo s.t. I can send in a figure and retrieve the logo fig.

import contextlib
from collections import defaultdict
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import re
import shutil
import subprocess
import sys
from urllib.request import urlopen

import config
from utils import seq_logo

def create_internal_dir(store_dir, replace_existing):
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

# todo: merge meme_out and mast_out into a single folder...?

def get_pdb_list(seq_cid_map, pdb_folder,
                 type='prosite_extract',
                 store_dir=config.store_dir,
                 replace_existing=False,
                 delete_intermediate_store=False,
                 ref_meme_txt=None, disable_tests=False):
    _dir = create_internal_dir(store_dir, replace_existing)
    if not disable_tests:
        _test_seq_cid_map(seq_cid_map)
    download_pdb_files(seq_cid_map, pdb_folder, replace_existing=False)
    if not disable_tests:
        _test_pdb_files_present(seq_cid_map, pdb_folder)
    seqs = extract_sequence(pdb_folder, seq_cid_map)
    write_to_file(seqs, _dir['fasta_fpath'])
    delete_short_seqs(_dir['fasta_fpath'], threshold=30)
    if not disable_tests:
        _test_fasta_match_pdb(_dir['fasta_fpath'], pdb_folder, seq_cid_map)

    if type == 'prosite_extract':
        run_meme(_dir['meme_out'], _dir['fasta_fpath'], replace_existing=False)
        if not disable_tests:
            _test_successful_meme(_dir['meme_out'])
        meme_txt_path = os.path.join(_dir['meme_out'], 'meme.txt')
        seq_motif_map = get_motif_diagram_meme(meme_txt_path)
    elif type == 'datalon':
        assert isinstance(ref_meme_txt, str)
        assert os.path.isfile(ref_meme_txt)
        run_mast(_dir['mast_out'], ref_meme_txt, _dir['fasta_fpath'],
                 replace_existing=True)
        if not disable_tests:
            _test_successful_mast(_dir['mast_out'])
        mast_txt_path = os.path.join(_dir['mast_out'], 'mast.txt')
        seq_motif_map = get_motif_diagram_mast(mast_txt_path)
    else:
        raise Exception
    seq_motif_map = adjust_motif_diagram(seq_motif_map)
    seq_motif_map = delete_gapped_motifs(seq_motif_map, _dir['fasta_fpath'])
    ptr_props = merge_ptr_props(seq_motif_map, seq_cid_map)

    if delete_intermediate_store:
        shutil.rmtree(store_dir)
    return ptr_props

def merge_ptr_props(seq_motif_map, pname_cid_map):
    ptr_props = defaultdict(dict)
    for pname, motif_pos in seq_motif_map.items():
        ptr_props[pname]['sno_markers'] = motif_pos
        ptr_props[pname]['cid'] = pname_cid_map[pname]
    return ptr_props

# ----------------------------------
# This downloads the .pdb files listed in pdb_list, from rcsb server.

def download_pdb_files(seq_cid_map, output_folder, replace_existing=True,
                           file_suffix='.pdb',
                           url_template='https://files.rcsb.org/view/{}.pdb'):
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    stored_pdb_files = set(os.listdir(output_folder))
    for pname in seq_cid_map.keys():
        pname = pname.lower()
        url = url_template.format(pname.strip())
        output_path = os.path.join(output_folder, pname+file_suffix)
        if replace_existing or pname+file_suffix not in stored_pdb_files:
            with contextlib.closing(urlopen(url)) as contents:
                with open(output_path, 'w') as output_file:
                    for line in contents:
                        output_file.write(line.decode("utf-8"))
    return True

# ----------------------------------
# For each pname, using the pname_cid map, we obtain the seq from the
# relevant cid, from the .pdb files. This is then written into a .fasta file,
# for processing using meme. We then delete sequences with length shorter
# than our 30-len analysis pattern.

def extract_sequence(pdb_folder, pname_cid_map, AA3_to_AA1=config.AA3_to_AA1):
    seqs = []
    for pname in pname_cid_map.keys():
        pname = pname.split(".")[0]
        cid = pname_cid_map[pname]
        assert len(cid) == 1
        seq = [">{}\n".format(pname)]
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
                if sno > current_atom_count:
                    num_spacer_res = sno - current_atom_count
                    seq.append("X" * num_spacer_res)
                    res = line[17:20].strip()
                    AA_single_letter = AA3_to_AA1[res]
                    seq.append(f"{AA_single_letter}")
                    current_atom_count = sno + 1
                elif sno == current_atom_count:
                    res = line[17:20].strip()
                    AA_single_letter = AA3_to_AA1[res]
                    seq.append(f"{AA_single_letter}")
                    current_atom_count += 1
            seq.append("\n")
        seq = "".join(seq)
        seqs.append(seq)
    return seqs

def delete_short_seqs(fasta_fname, threshold):
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
                output_str.append(tmp)
                output_str.append(line)
                continue
    with open(fasta_fname, 'w') as file:
        for line in output_str:
            file.write(line)
    return

def write_to_file(data, filename):
    assert isinstance(data, list)
    with open(filename, 'w') as file:
        file.writelines(data)
    return

# ----------------------------------
# meme/mast portion

def run_meme(meme_output, fasta_filename, replace_existing=False):
    # parallelise this eventually, with a merge_meme_output function.
    if replace_existing or not os.path.isdir(meme_output):
        command = f"meme -w 13 -protein -nmotifs 1 -mod anr -oc {meme_output} " \
            f"{fasta_filename}"
        subprocess.run(command, shell=True)
    return True

def run_mast(mast_output, meme_txt, fasta_filename, replace_existing=False):
    if replace_existing or not os.path.isdir(mast_output):
        command = f"mast -oc {mast_output} -mt 0.0001 {meme_txt} " \
            f"{fasta_filename}"
        subprocess.run(command, shell=True)
    return True

# ----------------------------------
# We obtain from meme output (meme.txt) the motif diagram, showing the
# location of the matched motifs for each pname sequence. We then adjust it
# such that the location are absolute, relative to the start of the sequence
# rather than from the end of the previous motif. Motifs with gaps in between
# are also deleted, mainly because the subsequent descriptor code assumes a
# continuous sequence for the len-30 analysis segment.

def get_motif_diagram_meme(meme_txt):
    seq_motif_map = dict()
    with open(meme_txt, 'r') as rfile:
        correct_segment = False
        in_area = False
        for line in rfile:
            if re.search("MEME-1 block diagrams", line):
                correct_segment = True
                continue
            if correct_segment and line.startswith("-------------   "):
                in_area = True
                continue
            if in_area and line.startswith("---------------------------"):
                break
            if in_area:
                seq_name, motif_diagram = line[:4], line[43:]
                raw_motif_positions = motif_diagram.split("_[1]_")
                motif_positions = []
                for pos in raw_motif_positions[:-1]:
                    if pos == "":
                        motif_positions.append(0)
                    else:
                        motif_positions.append(int(pos))
                seq_motif_map[seq_name] = motif_positions
    return seq_motif_map


def get_motif_diagram_mast(mast_txt):
    seq_motif_map = dict()
    count = 0
    with open(mast_txt, 'r') as rfile:
        correct_segment = False
        in_area = False
        for line in rfile:
            if line.startswith("SECTION II: MOTIF DIAGRAMS"):
                correct_segment = True
                continue
            if correct_segment and line.startswith("-------------       "):
                in_area = True
                continue
            if in_area and not line.strip():
                break
            if in_area:
                pname, motif_diagram = line[:4], line[45:]
                if "-[1]-" not in motif_diagram:
                    continue
                raw_motif_positions = motif_diagram.split("-[1]-")
                motif_positions = []
                count += len(raw_motif_positions) - 1
                for pos in raw_motif_positions[:-1]:
                    if pos == "":
                        motif_positions.append(0)
                    else:
                        motif_positions.append(int(pos))
                seq_motif_map[pname] = motif_positions
    return seq_motif_map

def adjust_motif_diagram(prev_map):
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

def delete_gapped_motifs(prev_map, fasta_fname):
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

# ----------------------------------
# This debugging function displays the sequence logo of the identified motifs.

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
            #     Might have been dropped because of gapped motif, etc
                pname = None
                continue
            for pos in motif_pos:
                if len(line) <= pos+13:
                    logging.error(
                        f"Fasta seq is shorter than pos+13, for pos in "
                        f"motif_pos. Fasta_seq: <{line}>, "
                        f"motif_pos: <{motif_pos}>, illegal pos: <{pos}>.")
                    raise Exception
                labelled_motifs.append(line[pos-1:pos+12])
    return labelled_motifs

def plot_labelled_motifs(seq_motif_map, fasta_filename):
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
    return

def _test_fasta_match_pdb(fasta_fname, pdb_folder, seq_cid_map,
                          AA3_to_AA1=config.AA3_to_AA1):
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
                            # print(pname)
                            # print(f_line)
                            # print(p_line)
                            # print(AA3_to_AA1[res])
                            # print(f_line[sno - 1])
                            # print("")
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
        pname, _remain = split_fname
        fnames.append(pname)
        # assert pname in seq_cid_map
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