import os
import re
import subprocess

from config import paths
from utils import generic


from tests.src import paths_test

def func():
    debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
    input_seqfile = paths_test.MEME_TEST_SEQ
    output = os.path.join(debug_folder, "output_meme_files")
    run_meme(input_seqfile, 30, output, num_p=7)
    meme_txt = os.path.join(output, 'meme.txt')
    from biopython_adapted import bio_interface
    bio_interface.parse_meme_file(meme_txt)

    # with open(meme_txt, 'r') as handle:
    #     motifs = meme.read(handle)
    # motif = motifs[0]meme_pathmeme_path
    # print(motif.length)
    # for i, j in motif.counts.items():
    #     print(i)
    #     print(j)
    # print(motif.counts.values())
    # print(type(motif.counts))


def run_meme(fasta_filename, motif_len, meme_output, num_p=1,
             meme_exec=paths.MEME_EXEC):
    assert motif_len >= 1
    assert isinstance(motif_len, int)
    generic.quit_if_missing(fasta_filename)
    generic.warn_if_exist(meme_output, filetype='folder')
    command = f"{meme_exec} -w {motif_len} -p {num_p} -protein -nmotifs 1 " \
        f"-mod anr -oc {meme_output} {fasta_filename}"
    subprocess.run(command, shell=True)
    _test_successful_meme(meme_output)
    return True

def _test_successful_meme(meme_out):
    assert os.path.isdir(meme_out)
    _meme_txt_path = os.path.join(meme_out, 'meme.txt')
    assert os.path.isfile(_meme_txt_path)
    return True

def run_mast(meme_txt, fasta_filename, mast_output, mast_exec=paths.MAST_EXEC):
    generic.quit_if_missing(meme_txt)
    generic.quit_if_missing(fasta_filename)
    generic.warn_if_exist(mast_output, filetype='folder')
    command = f"{mast_exec} -oc {mast_output} -mt 0.0001 {meme_txt} " \
        f"{fasta_filename}"
    subprocess.run(command, shell=True)
    _test_successful_mast(mast_output)
    return True

def extract_motifs_meme(input_txt, motif_len):
    pname_motif_raw = _get_motif_diagram(input_txt, source='meme')
    pname_motif_map = _adjust_motif_diagram(pname_motif_raw, motif_len)
    return pname_motif_map

def extract_motifs_mast(input_txt, motif_len):
    pname_motif_raw = _get_motif_diagram(input_txt, source='mast')
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


def _test_successful_mast(mast_out):
    assert os.path.isdir(mast_out)
    _mast_txt_path = os.path.join(mast_out, 'mast.txt')
    assert os.path.isfile(_mast_txt_path)


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
    pname_motif_map = dict()
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


# func()
# if __name__ == "__main__":
#     meme_txt = os.path.join(paths.ROOT, "meme.txt")
#
#     # fasta_filename = paths.TEMPLATE_SEQFILE
#     fasta_filename = os.path.join(paths.ROOT, "mg_100_filtered.fasta")
#     mast_output = os.path.join(paths.ROOT, "MY_OUTPUT")
#     run_mast(meme_txt, fasta_filename, mast_output, mast_exec=paths.MAST_EXEC)