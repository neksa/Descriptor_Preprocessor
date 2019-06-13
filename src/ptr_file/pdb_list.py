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
import logging
import pickle
import os
import re
from urllib.request import urlopen

from global_config import input_dir
import ptr_utils

pdb_folder = os.path.join(input_dir, "pdb_files")
fasta_filename = os.path.join(input_dir, "tmp_prosite_seqs.fasta")
prosite_extract_path = os.path.join(input_dir, "prosite_extract.txt")
meme_output = os.path.join(input_dir, "meme_out")

def get_pdb_list(pdb_folder, extract_file, pdb_list, output_path,
                 type='prosite_extract'):
    pname_cid_map = parse_prosite_extract(extract_file, pdb_list)
    download_pdb_files_for(pdb_list, pdb_folder, replace_existing=False)
    seqs = extract_sequence(pdb_folder, pname_cid_map)
    write_to_file(seqs, fasta_filename)
    delete_short_seqs(fasta_filename, threshold=30)
    # run meme
    command = f"meme -w 13 -protein -nmotifs 1 -mod anr -oc {meme_output} " \
        f"{fasta_filename}"
    meme_txt_path = os.path.join(meme_output, "meme.txt")
    seq_motif_map = get_motif_diagram(meme_txt_path)
    seq_motif_map = adjust_motif_diagram(seq_motif_map)
    seq_motif_map = delete_gapped_motifs(seq_motif_map, fasta_filename)
    with open(output_path, 'wb') as file:
        pickle.dump(seq_motif_map, file, -1)
    return True

# ----------------------------------
# We take as input a prosite extract, and obtain from it the mapping between
# pname and cid to be used.

def parse_prosite_extract(extract_file, pdb_list):
    with open(extract_file, 'r') as file:
        cids = []
        for line in file:
            if re.search("MOL_ID: 1", line):
                start, remain = line.split("MOL_ID: 1")
                value = re.search("CHAIN: ([A-Z1-9])", remain)
                if value is None:
                    logging.warning(f"Valid cid not found in prosite_extract "
                                    f"for line <{line}>.")
                cids.append(value.group(1))
    pname_cid_map = dict()
    for pname, cid in zip(pdb_list, cids):
        pname_cid_map[pname] = cid
    return pname_cid_map

# ----------------------------------
# This downloads the .pdb files listed in pdb_list, from rcsb server.
def download_pdb_files_for(pdb_list, output_folder, replace_existing=True,
                           file_suffix='.pdb',
                           url_template='https://files.rcsb.org/view/{}.pdb'):
    for pdb_code in pdb_list:
        url = url_template.format(pdb_code.strip())
        output_path = os.path.join(output_folder, pdb_code+file_suffix)
        if replace_existing or not os.path.isfile(output_folder):
            with urlopen(url) as contents:
                with open(output_path, 'w') as output_file:
                    for line in contents:
                        output_file.write(line.decode("utf-8"))
    return True

# ----------------------------------
# For each pname, using the pname_cid map, we obtain the seq from the
# relevant cid, from the .pdb files. This is then written into a .fasta file,
# for processing using meme. We then delete sequences with length shorter
# than our 30-len analysis pattern.

def extract_sequence(pdb_folder, pname_cid_map):
    seqs = []
    for pdb_filename in os.listdir(pdb_folder):
        pname = pdb_filename.split(".")[0]
        cid = pname_cid_map[pname]
        assert len(cid) == 1
        seq = [">{}\n".format(pname)]
        pdb_filepath = os.path.join(pdb_folder, pdb_filename)
        with open(pdb_filepath, 'r') as pdb_file:
            current_atom_count = 1
            for line in pdb_file:
                if line.startswith("MODEL        2"):
                    break
                if not line.startswith("ATOM"):
                    continue
                if line[21] != cid:
                    continue
                sno = int(line[22:26].strip())
                if sno > current_atom_count:
                    num_spacer_res = sno - current_atom_count
                    seq.append("X" * num_spacer_res)
                    current_atom_count = sno + 1
                elif sno == current_atom_count:
                    res = line[17:20].strip()
                    AA_single_letter = ptr_utils.AA_map[res]
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
# We obtain from meme output (meme.txt) the motif diagram, showing the
# location of the matched motifs for each pname sequence. We then adjust it
# such that the location are absolute, relative to the start of the sequence
# rather than from the end of the previous motif. Motifs with gaps in between
# are also deleted, mainly because the subsequent descriptor code assumes a
# continuous sequence for the len-30 analysis segment.

def get_motif_diagram(meme_txt):
    seq_motif_map = dict()
    count = 0
    with open(meme_txt, 'r') as rfile:
        correct_segment = False
        in_area = False
        for line in rfile:
            if line.startswith("	Motif DKDGBGTIDFEEF MEME-1 block diagrams"):
                correct_segment = True
                continue
            if correct_segment and line.startswith("-------------   "):
                in_area = True
                continue
            if in_area and line.startswith("---------------------------"):
                break
            if in_area:
                seq_name, motif_diagram = line[:8], line[43:]
                raw_motif_positions = motif_diagram.split("_[1]_")
                motif_positions = []
                count += len(raw_motif_positions) - 1
                for pos in raw_motif_positions[:-1]:
                    if pos == "":
                        motif_positions.append(0)
                    else:
                        motif_positions.append(int(pos))
                seq_motif_map[seq_name] = motif_positions
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

def _get_labelled_motifs(seq_motif_map, pname_cid_map, pdb_folder):
    labelled_motifs = []
    for seq_name, motif_pos in seq_motif_map.items():
        pdb_filename = os.path.join(pdb_folder, seq_name)
        cid = pname_cid_map[seq_name]
        sno_res_map = dict()

        with open(pdb_filename, 'r') as pdb_file:
            for line in pdb_file:
                if not line.startswith("ATOM"):
                    continue
                if line[21] != cid:
                    continue
                sno = int(line[22:26].strip())
                res = line[17:20]
                sno_res_map[sno] = ptr_utils.AA_map[res]

        for pos in motif_pos:
            motif = ""
            for sno in range(pos, pos+13):
                res = sno_res_map[sno] if sno in sno_res_map else "P"
                motif += res
            labelled_motifs.append(motif)
    return labelled_motifs