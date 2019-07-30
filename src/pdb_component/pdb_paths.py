import os

ROOT = os.path.dirname(__file__)
PDB_FILES = os.path.join(ROOT, "pdb_files")
PDB_PARSED = os.path.join(ROOT, "pdb_files_parsed")

PDB_FILES_SET = set(os.listdir(PDB_FILES))
PDB_PARSED_SET = set(os.listdir(PDB_PARSED))

HB_EXEC = os.path.join(ROOT, "parsers", "hb", "hb_calculator")

