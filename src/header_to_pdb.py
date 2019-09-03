"""
Takes in a header+m file, extract out the ACC_ID, match to .pdb id, download
.pdb files...? Call descr portion afterwards maybe.

"""
import find_cid_from_pname
from utils import logs

from config import paths
import preprocess
import os
from collections import Counter, OrderedDict

# todo: 1. From this: https://prosite.expasy.org/PS00164#TP, Get the aligned
#  sequences (the short pattern) in .fasta file, run meme on it to construct
#  a meme.txt, use mast on the meme txt, cut out the left-and-right portion
#  to make a length-30, and use that as my motif_pos.

# todo: meme takes forever. An alternative is to convert the aligned seqs
#  into a matrix, then call matrix2meme.

# todo: need to build a composition file using the full seq first.

import numpy as np
import pickle

from utils import generic, build_composition, build_meme_from_aligned, \
    get_pname_seq
import sys
# First, we download the prosite files, and derive our seed output_matrix:

file_folder = os.path.join(paths.ROOT, "Reports", "19082019", "PDOC00325")


full_seq_path = os.path.join(file_folder, "full_seq.fasta")
aligned_seq_path = os.path.join(file_folder, "aligner.txt")
output_matrix_path = os.path.join(file_folder, "output_matrix.txt")  #
# preprocess.run_prosite_aligned_cropped_flat(full_seq_path, aligned_seq_path,
#                                        output_matrix_path)
# preprocess.run_prosite_aligned_cropped(full_seq_path, aligned_seq_path,
#                                             output_matrix_path)
# import sys
# sys.exit()

# Then, we run converge_encoder, to encode 1. the input seed_matrix, 2. the
# prosite sequence file

# Then, we run on conv_rewrite, to obtain the output.1.matrix,
# composition.txt, and headers_m.txt

# Then, we take output.1.matrix, and attempt to create a logo using it. For
# that, we need to make composition.txt, which should be based on the
# sequences used to make the profile. However, since most matches will come
# from the aligned_seqs, it's probably k to just use those sequences to make
# our composition.

from converge.conv_interface import convert_conv_to_meme_full_num

converge_output_path = os.path.join(file_folder, "output.1.matrix")
composition_path = os.path.join(file_folder, "composition.txt")
output_meme_path = os.path.join(file_folder, "meme.txt")
# convert_conv_to_meme_full_num(converge_output_path, composition_path,
#                               output_meme_path)

# logo_path = os.path.join(file_folder, "logo.eps")
# import subprocess
#
# ceqlogo_exec = os.path.join(paths.ROOT, "src", "meme_suite", "meme", "ceqlogo")
# command = f"{ceqlogo_exec} -i1 {output_meme_path} -o {logo_path} -f EPS"
# subprocess.run(command, shell=True)



# After logo is done and checked, we then need to write a script that takes
# headers_m.txt, find the relevant .pdb id, download the .pdb files, align to
# find the chainid, then apply m to find the motif itself. Then, we build a
# descr based on that.
from utils import logs
logs.set_logging_level()


HEADERS_M_FILENAME = os.path.join(file_folder, "headers_m.txt")
acc_pos_map = dict()
with open(HEADERS_M_FILENAME, 'r') as file:
    for line in file:
        if line.startswith(">"):
            full_header, motif_pos = line.rsplit(",", maxsplit=1)
            acc_id = full_header.split("|")[1]
            acc_pos_map[acc_id] = int(motif_pos)

from utils import uniprot_id_converter


TMP_SEQ_FILE = os.path.join(file_folder, "tmp_fasta.txt")

from utils import download_fasta_given_pdb

# download_fasta_given_pdb.download_no_convert(list(acc_pos_map.keys()),
#                                              TMP_SEQ_FILE)

# tmp_composition_path = os.path.join(file_folder, "compos_align.txt")
# aligned_meme_path = os.path.join(file_folder, "aligned_meme.txt")
# build_composition.build(TMP_SEQ_FILE, tmp_composition_path)
# build_meme_from_aligned.build(aligned_seq_path, aligned_meme_path,
#                               tmp_composition_path)
# aligned_logo_path = os.path.join(file_folder, "aligned_logo.eps")
# command = f"{ceqlogo_exec} -i1 {aligned_meme_path} -o {aligned_logo_path} -f " \
#           f"EPS"
# subprocess.run(command, shell=True)

acc_seq_map = get_pname_seq.parse_raw(TMP_SEQ_FILE)
# print(acc_seq_map)


acc_ids = list(acc_pos_map.keys())
acc_pdb_map = uniprot_id_converter.convert("ACC", "PDB_ID", acc_ids)

pdb_seq_map = dict()
# because pdb => acc mapping may not be 1-1, we retain the original maps
mapped_pdb_acc = dict()
for acc_id, seq in acc_seq_map.items():
    if acc_id in acc_pdb_map:
        pdb_id = acc_pdb_map[acc_id]
        mapped_pdb_acc[pdb_id] = acc_id
        pdb_seq_map[pdb_id] = seq

pdb_cid_map = find_cid_from_pname.find(pdb_seq_map)

# acc_pdb_cid_map = {mapped_pdb_acc[pdb]: (pdb, cid) for pdb, cid in
#                    pdb_cid_map.items()}
# cropped_acc_list = list(acc_pdb_cid_map.keys())

pdb_cid_path = os.path.join(file_folder, "tmp_pdb_cid_map.pkl")

with open(pdb_cid_path, 'wb') as file:
    pickle.dump(pdb_cid_map, file, -1)

TMP_SEQ_PDB_FILE = os.path.join(file_folder, "tmp_pdb_fasta.fasta")

from pdb_component import pdb_paths
preprocess.create_seq(pdb_cid_path, TMP_SEQ_PDB_FILE, pdb_paths.PDB_FILES)
preprocess.filter_seq_file(TMP_SEQ_PDB_FILE, threshold=51)

import preprocess
from meme_suite import meme_interface
import motif_finder
from collections import defaultdict

cropped_seq_file = os.path.join(file_folder, "cropped_seq_tmp.fasta")

# preprocess.keep_only_acc(cropped_acc_list, TMP_SEQ_FILE, cropped_seq_file)
meme_interface.run_mast(output_meme_path, TMP_SEQ_PDB_FILE, paths.MEME_MAST_FOLDER)
mast_txt_path = os.path.join(paths.MEME_MAST_FOLDER, 'mast.txt')
pdb_motif_map = meme_interface.extract_motifs_mast(mast_txt_path, 50)
print(pdb_motif_map)
pdb_motif_map = motif_finder._delete_gapped_motifs(pdb_motif_map,
                                                   TMP_SEQ_PDB_FILE)
print(pdb_motif_map)
pdb_motif_pos = defaultdict(dict)
for pdb_id, motif_pos in pdb_motif_map.items():
    cid = pdb_cid_map[pdb_id]
    pdb_motif_pos[pdb_id]['sno_markers'] = motif_pos
    pdb_motif_pos[pdb_id]['cid'] = cid
output = os.path.join(file_folder, "stuff.pkl")
with open(output, 'wb') as file:
    pickle.dump(pdb_motif_pos, file, -1)

