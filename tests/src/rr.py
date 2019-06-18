import os
import config

path = os.path.join(config.ROOT, 'data', 'input', 'pdb_files')

for item in os.listdir(path):
    fpath = os.path.join(path, item)
    fpath_new = os.path.join(path, item.lower())
    os.rename(fpath, fpath_new)