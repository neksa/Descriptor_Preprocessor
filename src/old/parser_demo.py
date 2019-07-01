"""
Read file in provided path, return as single DF, specified in pdb_parser.py.

To Use:
Call Loader, supplying filepath. This loads filelines but does not parse them
yet. Calling the returned instance with parse_with(parser_in_str) returns the
DF. parser_in_str can be ATOMParser, HETATMParser, MODRESParser, HbondParser
for .pdb files, or hbParser for .hb2 files. Multiple parse_with() can be
called with a single Loader instance, reducing IO load.

If multiple files are required, make multiple calls using os.listdir and feed
output filenames to create multiple Loader instances, with directory added
to the front of filename to create full filepath.
"""

import os

import config
from descr.parsers import loader
from utils import logs

logs.set_logging_level()

input_dir = os.path.join(config.input_dir, "pdb_files")

if __name__ == "__main__":
    files = os.listdir(input_dir)
    pdb_files = []
    hb_files = []
    for file in files:
        if file.endswith('.pdb'):
            pdb_files.append(file)

    ATOM_lines = []
    HETATM_lines = []
    MODRES_lines = []
    HB_lines = []

    hb_lines = []
    start = False
    for pdb_file in pdb_files:
        pdb_loaded = loader.Loader(os.path.join(input_dir, pdb_file))
        pdb_ATOM = pdb_loaded.parse_with('ATOMParser')
        pdb_HETATM = pdb_loaded.parse_with('HETATMParser')
        pdb_MODRES = pdb_loaded.parse_with('MODRESParser')
        pdb_HB = pdb_loaded.parse_with('HbondParser')

        ATOM_lines.append(pdb_ATOM)
        HETATM_lines.append(pdb_HETATM)
        MODRES_lines.append(pdb_MODRES)
        HB_lines.append(pdb_HB)


    for hb_file in hb_files:
        hb_loaded = loader.Loader(os.path.join(input_dir, hb_file))
        hb = hb_loaded.parse_with('hbParser')
        hb_lines.append(hb)

    # logging.info(type(hb_lines[0]))   # pd.DataFrame
    # logging.info(type(ATOM_lines[0]))  # pd.DataFrame
    # logging.info(len(HETATM_lines[0])) # Number of files in directory ending
    # logging.info(len(MODRES_lines))    # with .pdb

    pass