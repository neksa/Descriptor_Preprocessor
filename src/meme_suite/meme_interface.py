import logging
import os
import re
import subprocess

from config import paths
from utils import generic


def run_meme_single(fasta_filename, motif_len, meme_output, num_p=1,
                    meme_exec=paths.MEME_EXEC):
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    generic.quit_if_missing(fasta_filename)
    generic.warn_if_exist(meme_output, filetype='folder')
    command = f"{meme_exec} -w {motif_len} -p {num_p} -protein -nmotifs 1 " \
              f"-mod oops -oc {meme_output} {fasta_filename}"
    return_code = subprocess.run(command, shell=True).returncode
    if return_code != 0:
        logging.error("run_meme_single() failed.")
        logging.error(f"Command: <{command}>")
        raise Exception


def run_meme(fasta_filename, motif_len, meme_output, num_p=1,
             meme_exec=paths.MEME_EXEC):
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    generic.quit_if_missing(fasta_filename)
    generic.warn_if_exist(meme_output, filetype='folder')
    command = f"{meme_exec} -w {motif_len} -p {num_p} -protein -nmotifs 1 " \
        f"-mod anr -oc {meme_output} {fasta_filename}"
    return_code = subprocess.run(command, shell=True).returncode
    if return_code != 0:
        logging.error("run_meme() failed.")
        logging.error(f"Command: <{command}>")
        raise Exception
    return True


def run_mast(meme_txt, fasta_filename, mast_output, mast_exec=paths.MAST_EXEC):
    generic.quit_if_missing(meme_txt)
    generic.quit_if_missing(fasta_filename)
    generic.warn_if_exist(mast_output, filetype='folder')
    command = f"{mast_exec} -oc {mast_output} -mt 0.0001 {meme_txt} " \
        f"{fasta_filename}"
    return_code = subprocess.run(command, shell=True).returncode
    if return_code != 0:
        logging.error("run_mast() failed.")
        logging.error(f"Command: <{command}>")
        raise Exception


def extract_motifs_meme(input_txt, motif_len):
    pname_motif_raw = _get_motif_diagram_meme_pdb(input_txt)
    pname_motif_map = _adjust_motif_diagram(pname_motif_raw, motif_len)
    return pname_motif_map


def extract_motifs_mast(input_txt, motif_len):
    pname_motif_raw = _get_motif_diagram_mast_pdb(input_txt)
    pname_motif_map = _adjust_motif_diagram(pname_motif_raw, motif_len)
    return pname_motif_map


def extract_motifs_mast_uniprot(input_txt, motif_len):
    pname_motif_raw = _get_motif_diagram_mast_uniprot(input_txt)
    pname_motif_map = _adjust_motif_diagram(pname_motif_raw, motif_len)
    return pname_motif_map


def _adjust_motif_diagram(prev_map, motif_len):
    pname_motif_map = dict()
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
        pname_motif_map[pname] = abs_motif_pos
        assert len(relative_motif_pos) == len(abs_motif_pos)
    return pname_motif_map


def _get_motif_diagram_meme_pdb(input_txt):
    """
    We obtain from meme/mast output (meme.txt/mast.txt) the motif diagram,
    showing the location of the matched motifs for each pname sequence. We
    then adjust it such that the locations are absolute, relative to the
    start of the sequence rather than from the end of the previous motif.
    Motifs with gaps in between are also deleted, mainly because the
    subsequent descriptor code assumes a continuous sequence for the len-30
    analysis segment.
    """
    pname_motif_map = dict()
    motif_diagram_sep = "_[1]_"
    with open(input_txt, 'r') as rfile:
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
                pname, motif_diagram = line[:4], line[43:]
                if motif_diagram_sep not in motif_diagram:
                    continue
                # Remove [1]_5 if first motif is right at the front.
                if motif_diagram.startswith("["):
                    motif_diagram = motif_diagram[4:]
                raw_motif_positions = motif_diagram.split(motif_diagram_sep)
                motif_positions = []
                for pos in raw_motif_positions[:-1]:
                    if pos == "":
                        motif_positions.append(0)
                    else:
                        motif_positions.append(int(pos))
                pname_motif_map[pname] = motif_positions
    return pname_motif_map


def _get_motif_diagram_mast_pdb(input_txt):
    """
    We obtain from meme/mast output (meme.txt/mast.txt) the motif diagram,
    showing the location of the matched motifs for each pname sequence. We
    then adjust it such that the locations are absolute, relative to the
    start of the sequence rather than from the end of the previous motif.
    Motifs with gaps in between are also deleted, mainly because the
    subsequent descriptor code assumes a continuous sequence for the len-30
    analysis segment.
    """
    pname_motif_map = dict()
    motif_diagram_sep = "-[1]-"
    with open(input_txt, 'r') as rfile:
        correct_segment = False
        in_area = False
        for line in rfile:
            if line.startswith("SECTION II: MOTIF DIAGRAMS"):
                correct_segment = True
                continue
            if correct_segment and line.startswith("-------------   "):
                in_area = True
                continue
            if in_area and not line.strip():
                break
            if in_area:
                pname, motif_diagram = line[:4], line[43:]
                if motif_diagram_sep not in motif_diagram:
                    continue
                # Remove [1]_5 if first motif is right at the front.
                if motif_diagram.startswith("["):
                    motif_diagram = motif_diagram[4:]
                raw_motif_positions = motif_diagram.split(motif_diagram_sep)
                motif_positions = []
                for pos in raw_motif_positions[:-1]:
                    if pos == "":
                        motif_positions.append(0)
                    else:
                        motif_positions.append(int(pos))
                pname_motif_map[pname] = motif_positions
    return pname_motif_map


def _get_motif_diagram_mast_uniprot(input_txt):
    """
    We obtain from meme/mast output (meme.txt/mast.txt) the motif diagram,
    showing the location of the matched motifs for each pname sequence. We
    then adjust it such that the locations are absolute, relative to the
    start of the sequence rather than from the end of the previous motif.
    Motifs with gaps in between are also deleted, mainly because the
    subsequent descriptor code assumes a continuous sequence for the len-30
    analysis segment.
    """
    pname_motif_map = dict()
    motif_diagram_sep = "-[1]-"
    with open(input_txt, 'r') as rfile:
        correct_segment = False
        in_area = False
        for line in rfile:
            if line.startswith("SECTION II: MOTIF DIAGRAMS"):
                correct_segment = True
                continue
            if correct_segment and line.startswith("-------------   "):
                in_area = True
                continue
            if in_area and not line.strip():
                break
            if in_area:
                try:
                    pname = line.split("|")[1]
                    motif_diagram = line[45:].strip()
                    if motif_diagram_sep not in motif_diagram:
                        continue
                    # Remove [1]_5 if first motif is right at the front.
                    if motif_diagram.startswith("["):
                        motif_diagram = motif_diagram[4:]
                    raw_motif_positions = motif_diagram.split(motif_diagram_sep)
                    motif_positions = []
                    for pos in raw_motif_positions[:-1]:
                        if pos == "":
                            motif_positions.append(0)
                        else:
                            motif_positions.append(int(pos))
                    pname_motif_map[pname] = motif_positions
                except:
                    print(line)
                    continue
    return pname_motif_map

# func()
# if __name__ == "__main__":
#     meme_txt = os.path.join(paths.ROOT, "meme.txt")
#
#     # fasta_filename = paths.TEMPLATE_SEQFILE
#     fasta_filename = os.path.join(paths.ROOT, "mg_100_filtered.fasta")
#     mast_output = os.path.join(paths.ROOT, "MY_OUTPUT")
#     run_mast(meme_txt, fasta_filename, mast_output, mast_exec=paths.MAST_EXEC)