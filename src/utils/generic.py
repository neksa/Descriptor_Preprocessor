import os
import contextlib
from urllib import request
import logging
import shutil

def warn_if_exist(path, filetype='file', remove=True):
    assert filetype in ('file', 'folder')
    if filetype == 'file':
        if os.path.isfile(path):
            logging.warning(f"File in <{path}> exists. Replacing.")
            if remove:
                os.remove(path)
    else:
        if os.path.isdir(path):
            logging.warning(f"Folder in <{path}> exists. Replacing.")
            if remove:
                shutil.rmtree(path)

def quit_if_missing(path, filetype='file'):
    assert filetype in ('file', 'folder')
    if filetype == 'file':
        if not os.path.isfile(path):
            logging.error(f"File in <{path}> missing, exiting.")
            raise Exception
    else:
        if not os.path.isdir(path):
            logging.error(f"Folder in <{path}> missing, exiting.")
            raise Exception

_aa_index = [('ALA', 'A'),
             ('CYS', 'C'),
             ('ASP', 'D'),
             ('GLU', 'E'),
             ('PHE', 'F'),
             ('GLY', 'G'),
             ('HIS', 'H'),
             ('HSE', 'H'),
             ('HSD', 'H'),
             ('ILE', 'I'),
             ('LYS', 'K'),
             ('LEU', 'L'),
             ('MET', 'M'),
             ('MSE', 'M'),
             ('ASN', 'N'),
             ('PRO', 'P'),
             ('GLN', 'Q'),
             ('ARG', 'R'),
             ('SER', 'S'),
             ('THR', 'T'),
             ('VAL', 'V'),
             ('TRP', 'W'),
             ('TYR', 'Y')]

AA3_to_AA1 = dict(_aa_index)
