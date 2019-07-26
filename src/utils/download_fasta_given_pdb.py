import contextlib
from urllib import request

from utils import uniprot_id_converter, generic


def download(pdb_list, output):
    pdb_acc_map = uniprot_id_converter.convert("PDB_ID", "ACC", pdb_list)
    generic.warn_if_exist(output)
    with open(output, 'w') as file:
        acc_unique = set(pdb_acc_map.values())
        for acc in acc_unique:
            url = f"https://www.uniprot.org/uniprot/{acc}.fasta"
            try:
                with contextlib.closing(request.urlopen(url)) as contents:
                    output = contents.read().decode("utf-8")
                    file.write(output)
            except:
                continue

