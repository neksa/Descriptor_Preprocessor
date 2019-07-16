import os

ROOT = "/home/melvin/Desktop/Predictive_for_efhand"
DATA = os.path.join(ROOT, 'data')
USER = os.path.join(DATA, 'user')
USER_INPUT = os.path.join(USER, 'input')
USER_OUTPUT = os.path.join(USER, 'output')
INTERNAL = os.path.join(DATA, 'internal')
DEBUG = os.path.join(DATA, 'debug')

PROSITE_EXTRACT = os.path.join(USER_INPUT, 'prosite_extract.txt')
IONCOM_EXTRACT = os.path.join(USER_INPUT, 'ioncom_extract.txt')
IONCOM_BINDING_SITES = os.path.join(USER_INPUT, 'ioncom_binding_sites.txt')
REF_MEME_TXT = os.path.join(USER_INPUT, 'ref_meme.txt')

PNAME_CID = os.path.join(INTERNAL, 'pname_cid_map.pkl')
PDB_FOLDER = os.path.join(INTERNAL, 'pdb_files')
MOTIF_POS = os.path.join(INTERNAL, 'motif_pos.pkl')
TMP = os.path.join(INTERNAL, 'tmp')
MEME_MAST_FOLDER = os.path.join(TMP, 'meme_mast')

FULL_SEQS = os.path.join(INTERNAL, 'seqs.fasta')
CONV_SEED_SEQS = os.path.join(INTERNAL, 'seed_seqs.fasta')
CONV_MEME_FILE = os.path.join(INTERNAL, 'cov_meme.txt')
TEMPLATE_SEQFILE = os.path.join(INTERNAL, 'seq_template.fasta')

BASH_EXEC = "/bin/bash"

# CONVERGE_EXEC = os.path.join(EXTERNAL, "converge")
#
# INPUT = os.path.join(DATA, 'input')
# EXTERNAL = os.path.join(ROOT, 'input')
#
#
# STORE = os.path.join(DATA, 'store')
# OUTPUT = os.path.join(DATA, 'output')
#
# PROSITE_INPUT = os.path.join(INPUT, "prosite_extract.txt")
# IONCOM_INPUT = os.path.join(INPUT, "ioncom.txt")
#
# PNAME_CID = os.path.join(STORE, "pname_cid_map.pkl")

# fasta_fpath = os.path.join(store_dir, "prosite_seqs.fasta")

# mast_out = os.path.join(store_dir, "mast_out")
# mast_txt_path = os.path.join(mast_out, "mast.txt")
