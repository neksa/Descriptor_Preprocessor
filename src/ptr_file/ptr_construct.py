"""
First, we obtain the list of relevant pbd code from the link, under PDB =>
see more:
https://prosite.expasy.org/PS00018

This is represented in pdb_list

Next, we selected [Detailed view], and copy the html code for the table into
a .txt file, in prosite_extract.txt.

Then, we download the .pdb files, for the pdb_list, from rcsb,
in download_pdb_files_for()

We extract the cid from prosite_extract.txt, and map the p_name to the cid,
in create_pname_cid_map()

We extract out the sequence for the relevant chain for each .pdb file,
to use eventually in meme, in extract_seq()

We then run meme, via

We extract from meme.txt the motif diagram, in get_motif_diagram()

# Because the motif diagram from meme.txt does not distinguish between sno,
# and because sno in the .pdb files are not always continuous or start from 0,
# we perform a matching of the motif diagram separators, to the sno in .pdb

We convert the motif_diagram to actual sno numbers, in convert_motif_diagram()

Finally, as a verification step, we take the labelled motifs based on their
motif diagrams, and produce the sequence logo,
in test_visualise_labelled_motifs()

It might make sense for us to parse every one of the relevant .pdb files,
and then store it in a .pkl form.
"""

# todo: need to run meme again. Also, write tests for fasta => pdb match,
#  before sending it to meme, maybe.

# todo: need to change atom parser as well.

# todo: perhaps try writing a proper ATOMParser, and trying it out on the
#  various pdb first.

import urllib.request
from ptr_utils import *
import os
import re








def download_pdb_files_for(pdb_list, replace_existing=True):
    for pdb_code in pdb_list:
        url = url_template.format(pdb_code.strip())
        output_filepath = store_dir.format("pdb_files/" + pdb_code)
        if replace_existing and not os.path.isfile(output_filepath):
            with urllib.request.urlopen(url) as contents:
                with open(output_filepath, 'w') as \
                        output_file:
                    for line in contents:
                        output_file.write(line.decode("utf-8"))

def parse_prosite_extract(filepath):
    # filepath = "prosite_extract.txt"
    with open(filepath, 'r') as file:
        extracted = []
        for line in file:
            if re.search("MOL_ID: 1", line):
                start, remain = line.split("MOL_ID: 1")
                value = re.search("CHAIN: ([A-Z1-9])", remain)
                if value is None:
                    print(line)
                extracted.append(value.group(1))
    return extracted

def create_pname_cid_map(pdb_list, prosite_extract):
    mapping = dict()
    for pname, cid in zip(pdb_list, prosite_extract):
        mapping[pname+".pdb"] = cid
    return mapping

def extract_seq(pdb_folder, pname_cid_map, fasta_filename):
    with open(fasta_filename, 'w') as fasta_file:
        for pdb_filename in os.listdir(pdb_folder):
            cid = pname_cid_map[pdb_filename]
            with open(os.path.join(pdb_folder, pdb_filename), 'r') \
                    as pdb_file:
                current_atom_count = 1
                fasta_file.write(">{}\n".format(pdb_filename))
                for line in pdb_file:
                    if not line.startswith("ATOM"):
                        continue
                    if line[21] != cid:
                        continue
                    sno = int(line[22:26].strip())
                    if sno > current_atom_count:
                        for i in range(sno-current_atom_count+1):
                            fasta_file.write("X")
                        current_atom_count = sno + 1
                    elif sno == current_atom_count:
                        res = line[17:20].strip()
                        AA_single_letter = AA_map[res]
                        fasta_file.write(f"{AA_single_letter}")
                        current_atom_count += 1
                    else:
                        continue
                fasta_file.write("\n")
    return

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
                separators = motif_diagram.split("_[1]_")
                _separators = []
                count += len(separators) - 1
                for sep in separators[:-1]:
                    if sep != "":
                        _separators.append(int(sep))
                    else:
                        _separators.append(0)
                seq_motif_map[seq_name] = _separators
    return seq_motif_map

def convert_motif_diagram(seq_motif_map):
    new_seq_motif_map = dict()
    for seq_name, _separators in seq_motif_map.items():
        separators = []
        curr_count = 1
        for i, sep in enumerate(_separators):
            if i == 0:
                new_sep = sep + 1
                curr_count += sep
                curr_count += 13
            else:
                new_sep = curr_count + sep
                curr_count += sep
                curr_count += 13
            separators.append(new_sep)
        new_seq_motif_map[seq_name] = separators
        assert len(_separators) == len(separators)
    return new_seq_motif_map

def _test_visualise_labelled_motifs(seq_motif_map, pname_cid_map, pdb_folder):
    extracted_motifs = []
    for seq_name, separators in seq_motif_map.items():
        pdb_filename = os.path.join(pdb_folder, seq_name)
        cid = pname_cid_map[seq_name]
        sno_to_res_map = dict()

        with open(pdb_filename, 'r') as pdb_file:
            for line in pdb_file:
                if not line.startswith("ATOM"):
                    continue
                if line[21] != cid:
                    continue
                sno = int(line[22:26].strip())
                res = line[17:20]
                sno_to_res_map[sno] = AA_map[res]
        # print(pdb_filename)
        # print(sno_to_res_map)
        # print("")

        for sep in separators:
            motif = ""
            for i in range(sep, sep+13):
                if i in sno_to_res_map:
                    motif += sno_to_res_map[i]
                else:
                    motif += "P"
            extracted_motifs.append(motif)
    return extracted_motifs

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

def delete_gapped_motifs(seq_motif_map, fasta_fname):
    new_seq_motif_map = dict()
    selected_pname = None
    with open(fasta_fname, 'r') as file:
        for line in file:
            if line.startswith(">"):
                selected_pname = line[1:].strip()
                continue
            if selected_pname and selected_pname in seq_motif_map:
                cropped_separators = []
                orig_separators = seq_motif_map[selected_pname]
                for sep in orig_separators:
                    try:
                        motif = line[sep-8:sep+22]
                        for char in motif:
                            if char == "X":
                                break
                        else:
                            cropped_separators.append(sep)
                    except:
                        print("delete_gapped_motifs")
                        print(selected_pname)
                        print("")
                        continue
                if cropped_separators:
                    new_seq_motif_map[selected_pname] = cropped_separators
    return new_seq_motif_map


from global_config import input_dir

# todo: ignore all motifs that isn't continuous in the sno.

prosite_extract_path = os.path.join(input_dir, "prosite_extract.txt")
pdb_folder = os.path.join(input_dir, "pdb_files")
fasta_filename = os.path.join(input_dir, "tmp_prosite_seqs.fasta")

prosite_extract = parse_prosite_extract(prosite_extract_path)
prosite_extract = []
pname_cid_map = create_pname_cid_map(pdb_list, prosite_extract)
#
# extract_seq(pdb_folder, pname_cid_map, fasta_filename)
# delete_short_seqs(fasta_filename, threshold=30)

# import subprocess
# command = f"meme -w 13 -protein -nmotifs 1 -mod anr -p 7 {fasta_filename}"
# subprocess.run(command, shell=True)

# import sys
# sys.exit()

seq_motif_map = get_motif_diagram(os.path.join(input_dir, "meme.txt"))
seq_motif_map = convert_motif_diagram(seq_motif_map)
seq_motif_map = delete_gapped_motifs(seq_motif_map, fasta_filename)
print(seq_motif_map)

assembed_ptr_data = dict()
for pname, separators in seq_motif_map.items():
    to_add = dict()
    to_add['sno_markers'] = separators
    to_add['cid'] = pname_cid_map[pname]
    assembed_ptr_data[pname] = to_add


#
# extracted_motifs = _test_visualise_labelled_motifs(seq_motif_map,
#                                                  pname_cid_map, pdb_folder)
# print(extracted_motifs)
import matplotlib.pyplot as plt
from utils.seq_logo import Logo
import pickle
import numpy as np

print(os.listdir("./"))

with open("assembed_ptr_data.pkl", 'wb') as file:
    pickle.dump(assembed_ptr_data, file, -1)

# with open("test.pkl", 'rb') as file:
#     extracted_motifs = pickle.load(file)
# print(extracted_motifs)
import sys
sys.exit()
plt.figure()

converted_motifs = []
extracted_motifs = np.array([list(i) for i in extracted_motifs])
# for i in extracted_motifs:
#     print(len(i))
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

Logo(converted_motifs, -1, convert_AA3=False)
plt.show()
#
#
#










# def match_motif_diagram_to_pdb_sno(pdb_folder, pname_cid_map, seq_motif_map):
#     # replace fasta_outdir with fasta_filename
#     new_seq_motif_map = dict()
#     for pdb_filename in os.listdir(pdb_folder):
#         cid = pname_cid_map[pdb_filename]
#         with open(os.path.join(pdb_folder, pdb_filename), 'r') \
#                 as pdb_file:
#             current_atom_count = 1
#             fasta_file.write(">{}\n".format(pdb_filename))
#             for line in pdb_file:
#                 if not line.startswith("ATOM"):
#                     continue
#                 if line[21] != cid:
#                     continue
#                 sno = int(line[22:26].strip())
#                 if sno > current_atom_count:
#                     for i in range(sno-current_atom_count):
#                         fasta_file.write("X")
#                     current_atom_count = sno + 1
#                 elif sno == current_atom_count:
#                     res = line[17:20].strip()
#                     AA_single_letter = AA_map[res]
#                     fasta_file.write(f"{AA_single_letter}")
#                     current_atom_count += 1
#                 else:
#                     continue
#     return