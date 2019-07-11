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

def find(pname_cid_map, process,
         num_p, ref_meme_txt, mast_meme_folder, seq_file):
    assert process in ('meme', 'mast')
    if process == 'mast':
        assert isinstance(ref_meme_txt, str)
        assert os.path.isfile(ref_meme_txt)
    seq_motif_map = _build_seq_motif_map(process,
                                         mast_meme_folder,
                                         seq_file, num_p,
                                         ref_meme_txt=ref_meme_txt)
    seq_motif_map = _adjust_motif_diagram(seq_motif_map)
    seq_motif_map = _delete_gapped_motifs(seq_motif_map, seq_file)
    motif_pos = _assemble_motif_pos(seq_motif_map, pname_cid_map)
    return motif_pos

def _run_meme(meme_output, fasta_filename, num_p):
    generic.quit_if_missing(fasta_filename)
    generic.warn_if_exist(meme_output, filetype='folder', remove=True)
    command = f"meme -w 13 -p {num_p} -protein -nmotifs 1 -mod anr " \
        f"-oc {meme_output} {fasta_filename}"
    subprocess.run(command, shell=True)
    _test_successful_meme(meme_output)
    return True

def _run_mast(mast_output, meme_txt, fasta_filename):
    generic.quit_if_missing(meme_txt)
    generic.quit_if_missing(fasta_filename)
    generic.warn_if_exist(mast_output, filetype='folder', remove=True)
    command = f"mast -oc {mast_output} -mt 0.0001 {meme_txt} " \
        f"{fasta_filename}"
    subprocess.run(command, shell=True)
    _test_successful_mast(mast_output)
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

def _build_seq_motif_map(process, tmp_output_folder, seq_file, num_p,
                         ref_meme_txt=None):
    if process == 'meme':
        _run_meme(tmp_output_folder, seq_file, num_p)
        _test_successful_meme(tmp_output_folder)
        meme_txt_path = os.path.join(tmp_output_folder, 'meme.txt')
        seq_motif_map = _get_motif_diagram(meme_txt_path, 'meme')
    elif process == 'mast':
        print(tmp_output_folder)
        print(ref_meme_txt)
        print(seq_file)
        print("\n")
        _run_mast(tmp_output_folder, ref_meme_txt, seq_file)
        _test_successful_mast(tmp_output_folder)
        mast_txt_path = os.path.join(tmp_output_folder, 'mast.txt')
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
