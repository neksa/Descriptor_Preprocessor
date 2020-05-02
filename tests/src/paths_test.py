import os

from config import paths

ROOT = paths.ROOT
ROOT_TEST = os.path.join(ROOT, 'tests')
DATA = os.path.join(ROOT_TEST, 'data')

USER_INPUT = os.path.join(DATA, 'input')
USER_OUTPUT = os.path.join(DATA, 'output')
INTERNAL = os.path.join(DATA, 'internal')
DEBUG = os.path.join(DATA, 'debug')
REF = os.path.join(DATA, 'ref')

MEME_SUITE = paths.MEME_SUITE
MEME_EXEC = paths.MEME_EXEC
MAST_EXEC = paths.MAST_EXEC
PDB_FOLDER = paths.PDB_FOLDER

PROSITE_EXTRACT = os.path.join(USER_INPUT, 'prosite_extract.txt')
IONCOM_EXTRACT = os.path.join(USER_INPUT, 'ioncom_extract.txt')
IONCOM_BINDING_SITES = os.path.join(USER_INPUT, 'ioncom_binding_sites.txt')
REF_MEME_TXT = os.path.join(USER_INPUT, 'ref_meme.txt')

TMP_FILE_TEMPLATE = os.path.join(DEBUG, 'output_{}.txt')

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
REF_RUN_PROSITE_MEME = os.path.join(REF, 'run_prosite_meme.pkl')
REF_RUN_PROSITE_MAST = os.path.join(REF, 'run_prosite_mast.pkl')
REF_RUN_IONCOM_MAST = os.path.join(REF, 'run_ioncom_meme.pkl')

# TestMemeConvConversion
ORIG_MEME_FOR_CONVERT = os.path.join(USER_INPUT, 'orig_meme.txt')
REF_CONV_MATRIX = os.path.join(REF, 'meme_to_conv_matrix.txt')
REF_CONV_COMPOSITION = os.path.join(REF, 'meme_to_conv_composition.txt')
REF_MEME_FROM_CONV = os.path.join(REF, 'converted_meme.txt')

# TestMemeSuite
MEME_TEST_SEQ = os.path.join(USER_INPUT, 'meme_test_seq.fasta')
MAST_TEST_DIAG = os.path.join(USER_INPUT, 'mast_test_diagram.txt')
REF_MEME_COUNTS = os.path.join(REF, 'meme_counts.pkl')
REF_MAST_DIAGRAMS = os.path.join(REF, 'mast_diagrams.pkl')

# TestConverge
CONV_TEST_SEQ = os.path.join(USER_INPUT, 'conv_test_seq.fasta')
INPUT_CONV_SEED_SEQS = os.path.join(USER_INPUT, 'input_conv_seed_seqs.fasta')
REF_CONV_WITH_BIND = os.path.join(REF, 'run_conv_bind.pkl')
REF_CONV_WITH_SEED = os.path.join(REF, 'run_conv_seed.pkl')

# TestGetPnameSeq
UNIPROT_SEQ = os.path.join(USER_INPUT, "uniprot_template_seq.fasta")
UNIREF_SEQ = os.path.join(USER_INPUT, "uniref_template_seq.fasta")
REF_UNIPROT_PNAME_SEQ = os.path.join(REF, 'uniprot_pname_seq.pkl')
REF_UNIREF_PNAME_SEQ = os.path.join(REF, 'uniref_pname_seq.pkl')



# PNAME_CID = os.path.join(INTERNAL, 'pname_cid_map.pkl')
# # PDB_FOLDER = os.path.join(ROOT, 'data', 'internal', 'pdb_files')
# MOTIF_POS = os.path.join(INTERNAL, 'motif_pos.pkl')

# MEME_MAST_FOLDER = os.path.join(TMP, 'meme_mast')
#
# FULL_SEQS = os.path.join(INTERNAL, 'seqs.fasta')
# CONV_SEED_SEQS = os.path.join(INTERNAL, 'seed_seqs.fasta')

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
