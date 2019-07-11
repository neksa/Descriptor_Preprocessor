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

def download_pdb_files(seq_cid_map,
                        output_folder,
                        replace_existing=True,
                        file_suffix='.pdb',
                        url_template='https://files.rcsb.org/view/{}.pdb'):
    # This downloads the .pdb files listed in pdb_list, from rcsb server.
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    stored_pdb_files = set(os.listdir(output_folder))
    for pname in seq_cid_map.keys():
        pname = pname.lower()
        url = url_template.format(pname.strip())
        output_path = os.path.join(output_folder, pname+file_suffix)
        if replace_existing or pname+file_suffix not in stored_pdb_files:
            with contextlib.closing(request.urlopen(url)) as contents:
                with open(output_path, 'w') as output_file:
                    output_file.write(contents.read().decode("utf-8"))
                    # for line in contents:
                    #     output_file.write(line.decode("utf-8"))

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
