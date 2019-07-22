import os

ROOT = "/home/melvin/Desktop/Predictive_for_efhand"
ROOT_TEST = os.path.join(ROOT, 'tests')
DATA = os.path.join(ROOT_TEST, 'data')
USER = os.path.join(DATA, 'user')
USER_INPUT = os.path.join(USER, 'input')
USER_OUTPUT = os.path.join(USER, 'output')
INTERNAL = os.path.join(DATA, 'internal')
REF = os.path.join(DATA, 'ref')

TMP = os.path.join(INTERNAL, 'tmp')
TMP_FILE = os.path.join(TMP, 'tmp_file')
TMP_FILE_TEMPLATE = os.path.join(TMP, 'tmp_file_{}')

PROSITE_EXTRACT = os.path.join(USER_INPUT, 'prosite_extract.txt')
IONCOM_EXTRACT = os.path.join(USER_INPUT, 'ioncom_extract.txt')
IONCOM_BINDING_SITES = os.path.join(USER_INPUT, 'ioncom_binding_sites.txt')
REF_MEME_TXT = os.path.join(USER_INPUT, 'ref_meme.txt')

# TestParseExtract
REF_PROSITE_EXTRACT_PNAME_CID = os.path.join(REF, 'prosite_extract.pkl')
REF_IONCOM_EXTRACT_PNAME_CID = os.path.join(REF, 'ioncom_extract.pkl')

# TestCreateSeq
REF_CREATE_SEQ_1 = os.path.join(REF, 'create_seq_1.fasta')
REF_CREATE_SEQ_2 = os.path.join(REF, 'create_seq_2.fasta')

# TestFindMotif
REF_FIND_MOTIF_1 = os.path.join(REF, 'find_motif_1.pkl')
REF_FIND_MOTIF_2 = os.path.join(REF, 'find_motif_2.pkl')

# TestRunAll
REF_RUN_ALL_1 = os.path.join(REF, 'run_all_1.pkl')
REF_RUN_ALL_2 = os.path.join(REF, 'run_all_2.pkl')
REF_RUN_ALL_3 = os.path.join(REF, 'run_all_3.pkl')

# TestMemeConvConversion
ORIG_MEME_FOR_CONVERT = os.path.join(USER_INPUT, 'orig_meme.txt')
REF_CONV_MATRIX = os.path.join(REF, 'meme_to_conv_matrix.txt')
REF_CONV_COMPOSITION = os.path.join(REF, 'meme_to_conv_compostion.txt')
REF_MEME_FROM_CONV = os.path.join(REF, 'converted_meme.txt')

PNAME_CID = os.path.join(INTERNAL, 'pname_cid_map.pkl')
PDB_FOLDER = os.path.join(ROOT, 'data', 'internal', 'pdb_files')
MOTIF_POS = os.path.join(INTERNAL, 'motif_pos.pkl')
TMP = os.path.join(INTERNAL, 'tmp')
MEME_MAST_FOLDER = os.path.join(TMP, 'meme_mast')

FULL_SEQS = os.path.join(INTERNAL, 'seqs.fasta')
CONV_SEED_SEQS = os.path.join(INTERNAL, 'seed_seqs.fasta')

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
