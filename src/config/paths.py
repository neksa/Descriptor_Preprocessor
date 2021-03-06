import os


ROOT = "/".join(os.path.dirname(__file__).split("/")[:-2])
BASH_EXEC = "/bin/bash"

# ROOT = "/home/yincp/Desktop/work/Descriptor_Preprocessor"
DATA = os.path.join(ROOT, 'data')
SRC = os.path.join(ROOT, 'src')
EXTERNAL = os.path.join(ROOT, 'external')
TRASH = os.path.join(ROOT, 'trash')

USER_INPUT = os.path.join(DATA, 'input')
USER_OUTPUT = os.path.join(DATA, 'output')
INTERNAL = os.path.join(DATA, 'internal')
DEBUG = os.path.join(DATA, 'debug')

TMP_FILE_TEMPLATE = os.path.join(DEBUG, 'tmp_{}')

MEME_SUITE = os.path.join(SRC, 'meme_suite')
MEME_EXEC = os.path.join(MEME_SUITE, 'meme', 'bin', 'meme')
MAST_EXEC = os.path.join(MEME_SUITE, 'meme', 'bin', 'mast')
FASTA_2_MARKOV_EXEC = os.path.join(MEME_SUITE, 'meme', 'src',
                                   'fasta-get-markov')
MATRIX_2_MEME_EXEC = os.path.join(MEME_SUITE, 'meme', 'scripts',
                                  'matrix2meme')
PDB_FOLDER = os.path.join(INTERNAL, 'pdb_files')

PROSITE_EXTRACT = os.path.join(USER_INPUT, 'prosite_extract.txt')
IONCOM_EXTRACT = os.path.join(USER_INPUT, 'ioncom_extract.txt')
IONCOM_BINDING_SITES = os.path.join(USER_INPUT, 'ioncom_binding_sites.txt')
REF_MEME_TXT = os.path.join(USER_INPUT, 'ref_meme.txt')

PNAME_CID = os.path.join(INTERNAL, 'pname_cid_map.pkl')
MOTIF_POS = os.path.join(INTERNAL, 'motif_pos.pkl')
MEME_MAST_FOLDER = os.path.join(INTERNAL, 'meme_mast')

FULL_SEQS = os.path.join(INTERNAL, 'seqs.fasta')
CONV_SEED_SEQS = os.path.join(INTERNAL, 'seed_seqs.fasta')
CONV_MEME_FILE = os.path.join(INTERNAL, 'cov_meme.txt')
TEMPLATE_SEQFILE = os.path.join(INTERNAL, 'seq_template.fasta')

PID_PDB_MAP = os.path.join(USER_OUTPUT, 'pid_pdb_map.pkl')

CONV_FOLDER = os.path.join(SRC, 'converge')
UNIPROT_BINARY = os.path.join(INTERNAL, 'uniprot_reviewed_all.binary')
RCSB_SEQS_FASTA = os.path.join(USER_INPUT, "rcsb_seqs_full.txt")
CONV_INPUT_MATRIX = os.path.join(CONV_FOLDER, 'input_matrix.txt')
CONV_OUTPUT = os.path.join(CONV_FOLDER, 'converged_matrix.txt')
CONV_COMPOSITION = os.path.join(CONV_FOLDER, 'composition.txt')
CONV_EXEC = os.path.join(CONV_FOLDER, 'calculator')

PROSITE_ALIGNED_SEQS = os.path.join(USER_INPUT,
                                    'prosite_enolase_aligned_seqs.txt')
PROSITE_ENOLASE_SEQS = os.path.join(USER_INPUT,
                                    'uniprot_prosite_enolase_seqs.fasta')

PRELOADED_PDB_FOLDER = os.path.join(INTERNAL, 'pdb_files_parsed')

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
