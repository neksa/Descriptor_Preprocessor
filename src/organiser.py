import os
import shutil

for filename in os.listdir("../All_PDB"):
    if not filename.endswith(".pdb"):
        continue

    a, b = filename[0], filename[1]
    if not os.path.isdir(f"../All_PDB/{a}"):
        os.mkdir(f"../All_PDB/{a}")
    if not os.path.isdir(f"../All_PDB/{a}/{b}"):
        os.mkdir(f"../All_PDB/{a}/{b}")
    shutil.move(f"../All_PDB/{filename}", f"../All_PDB/{a}/{b}/{filename}")
