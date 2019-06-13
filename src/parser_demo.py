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

# todo: spawn multiprocessing for this

import os
import logging

from global_config import input_dir
# This module should not be imported, otherwise .loader should be used.
# See intra-package referencing.
from descr.parsers.loader import Loader
from utils.logs import set_logging_level

set_logging_level()

input_dir = os.path.join(input_dir, "pdb_files")

if __name__ == "__main__":
    files = os.listdir(input_dir)
    pdb_files = []
    hb2_files = []
    for file in files:
        if file.endswith('.pdb'):
            pdb_files.append(file)
        # if file == '1B4C':
        #     pdb_files.append(file)
        # if file.endswith('.pdb'):
        #     if file == '1exr.pdb':
        #         pdb_files.append(file)
        # elif file.endswith('.hb2'):
        #     if file == '1exr.hb2':
        #         hb2_files.append(file)

    ATOM_lines = []
    HETATM_lines = []
    MODRES_lines = []
    HB_lines = []

    hb2_lines = []
    # 4BYF
    start = False
    # 6FEH check cid, and pdb
    # {512, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 513, 514, 532, 27, 28, 515, 29, 30, 31, 548, 516, 549, 32, 33, 34, 517, 35, 36, 37, 38, 518, 39, 40, 41, 42, 519, 44, 45, 46, 47, 520, 49, 50, 51, 52, 521, 54, 55, 56, 57, 522, 59, 60, 61, 62, 523, 64, 65, 66, 67, 524, 69, 70, 71, 72, 525, 74, 75, 76, 77, 526, 79, 80, 81, 82, 527, 84, 85, 86, 87, 528, 89, 90, 91, 92, 529, 94, 95, 96, 97, 530, 99, 100, 101, 102, 531, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 68, 103, 43, 48, 53, 58, 310, 311, 312, 313, 314, 315, 63, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 73, 392, 393, 78, 83, 88, 93, 98, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511}
    # 237

    for pdb_file in pdb_files:
        pdb_loaded = Loader(os.path.join(input_dir, pdb_file))
        pdb_ATOM = pdb_loaded.parse_with('ATOMParser')
        pdb_HETATM = pdb_loaded.parse_with('HETATMParser')
        pdb_MODRES = pdb_loaded.parse_with('MODRESParser')
        pdb_HB = pdb_loaded.parse_with('HbondParser')

        print(set(pdb_ATOM.sno))
        print(len(set(pdb_ATOM.sno)))
        print("")

        ATOM_lines.append(pdb_ATOM)
        HETATM_lines.append(pdb_HETATM)
        MODRES_lines.append(pdb_MODRES)
        HB_lines.append(pdb_HB)


    for hb2_file in hb2_files:
        hb2_loaded = Loader(os.path.join(input_dir, hb2_file))
        hb2 = hb2_loaded.parse_with('Hb2Parser')
        hb2_lines.append(hb2)

    # print(ATOM_lines)
    # print(set(ATOM_lines[0].sno))
    # print(len(set(ATOM_lines[0].sno)))
    # logging.info(HB_lines[0])
    # logging.info(hb2_lines[0])

    # logging.info(type(hb2_lines[0]))   # pd.DataFrame
    # logging.info(type(ATOM_lines[0]))  # pd.DataFrame
    # logging.info(len(HETATM_lines[0])) # Number of files in directory ending
    # logging.info(len(MODRES_lines))    # with .pdb

    pass