import os

src = os.path.dirname(__file__)
ROOT = src.rsplit(os.sep, 1)[0]
DATA = os.path.join(ROOT, 'data')
INPUT = os.path.join(DATA, 'input')
EXTERNAL = os.path.join(ROOT, 'input')

PDB_FOLDER = os.path.join(INPUT, 'pdb_files')
STORE = os.path.join(DATA, 'store')
OUTPUT = os.path.join(DATA, 'output')

PROSITE_INPUT = os.path.join(INPUT, "prosite_extract.txt")
IONCOM_INPUT = os.path.join(INPUT, "ioncom.txt")

PNAME_CID = os.path.join(STORE, "pname_cid_map.pkl")

# fasta_fpath = os.path.join(store_dir, "prosite_seqs.fasta")

# mast_out = os.path.join(store_dir, "mast_out")
# mast_txt_path = os.path.join(mast_out, "mast.txt")
CONVERGE_EXEC = os.path.join(EXTERNAL, "converge")