"""
Read file in provided path, return as single DF, specified in pdb_parser.py.

To Use:
Call Loader, supplying filepath. This loads filelines but does not parse them
yet. Calling the returned instance with parse_with(parser_in_str) returns the
DF. parser_in_str can be ATOMParser, HETATMParser, MODRESParser, HbondParser
for .pdb files, or Hb2Parser for .hb2 files. Multiple parse_with() can be
called with a single Loader instance, reducing IO load.

If multiple files are required, make multiple calls using os.listdir and feed
output filenames to create multiple Loader instances, with directory added
to the front of filename to create full filepath.
"""
# DEPRECATED

import os
import logging

from global_config import input_dir
# This module should not be imported, otherwise .loader should be used.
# See intra-package referencing.
from descr.parsers.loader import Loader
from utils.logs import set_logging_level

set_logging_level()

if __name__ == "__main__":
    files = os.listdir(input_dir)
    pdb_files = []
    hb2_files = []
    for file in files:
        if file.endswith('.pdb'):
            if file == '1exr.pdb':
                pdb_files.append(file)
        elif file.endswith('.hb2'):
            if file == '1exr.hb2':
                hb2_files.append(file)

    ATOM_lines = []
    HETATM_lines = []
    MODRES_lines = []
    HB_lines = []

    hb2_lines = []

    for pdb_file in pdb_files:
        print(pdb_file)
        pdb_loaded = Loader(os.path.join(input_dir, pdb_file))
        pdb_ATOM = pdb_loaded.parse_with('ATOMParser')
        pdb_HETATM = pdb_loaded.parse_with('HETATMParser')
        pdb_MODRES = pdb_loaded.parse_with('MODRESParser')
        pdb_HB = pdb_loaded.parse_with('HbondParser')

        ATOM_lines.append(pdb_ATOM)
        HETATM_lines.append(pdb_HETATM)
        MODRES_lines.append(pdb_MODRES)
        HB_lines.append(pdb_HB)

    for hb2_file in hb2_files:
        hb2_loaded = Loader(os.path.join(input_dir, hb2_file))
        hb2 = hb2_loaded.parse_with('Hb2Parser')
        hb2_lines.append(hb2)

    # logging.info(HB_lines[0])
    # logging.info(hb2_lines[0])

    # logging.info(type(hb2_lines[0]))   # pd.DataFrame
    # logging.info(type(ATOM_lines[0]))  # pd.DataFrame
    # logging.info(len(HETATM_lines[0])) # Number of files in directory ending
    # logging.info(len(MODRES_lines))    # with .pdb

    pass