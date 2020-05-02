
import os
import logging
import traceback

import loaders
from utils import logs
import pickle
from config import paths


def preload(pdb_folder, output_folder):
    assert os.path.isdir(pdb_folder)
    print(len(os.listdir(pdb_folder)))
    for i, filename in enumerate(sorted(os.listdir(pdb_folder))):
        pdb_id = filename.split(".")[0]
        print(f"{i}: {pdb_id}")
        filepath = os.path.join(pdb_folder, filename)
        try:
            file_data = loaders._load_data(filepath)
        except Exception as e:
            logging.error(f"_load_data() fails for file {filepath}. Skipping.")
            logging.error(f"Traceback: <{traceback.format_exc()}>")
            logging.error(f"Error_msg: <{e}>\n\n")
            continue
        output_path = os.path.join(output_folder, pdb_id+".pkl")
        with open(output_path, 'wb') as file:
            pickle.dump(file_data, file, -1)


def preload_update(pdb_folder, output_folder):
    assert os.path.isdir(pdb_folder)
    assert os.path.isdir(output_folder)
    print(len(os.listdir(pdb_folder)))
    output_filenames = set(os.listdir(output_folder))
    for i, filename in enumerate(sorted(os.listdir(pdb_folder))):
        pdb_id = filename.split(".")[0]
        if pdb_id+".pkl" in output_filenames:
            continue
        print(f"{i}: {filename}")
        filepath = os.path.join(pdb_folder, filename)
        try:
            file_data = loaders._load_data(filepath)
        except Exception as e:
            logging.error(f"_load_data() fails for file {filepath}. Skipping.")
            logging.error(f"Traceback: <{traceback.format_exc()}>")
            logging.error(f"Error_msg: <{e}>\n\n")
            continue
        output_path = os.path.join(output_folder, pdb_id + ".pkl")
        with open(output_path, 'wb') as file:
            pickle.dump(file_data, file, -1)

if __name__ == "__main__":
    logs.set_logging_level()
    preload(paths.PDB_FOLDER, paths.PRELOADED_PDB_FOLDER)