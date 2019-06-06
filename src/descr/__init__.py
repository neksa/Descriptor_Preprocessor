import os

current_dir = os.path.dirname(__file__)
new_dir = current_dir.rsplit(os.sep, 1)[0]
os.chdir(new_dir)