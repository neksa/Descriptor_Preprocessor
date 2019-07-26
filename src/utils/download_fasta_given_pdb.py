import contextlib
import os
from urllib import request

from utils import uniprot_id_converter, generic

# def query_genename(genenames):
#     id_pdb_map = _query(genenames, "GENENAME")
#     return id_pdb_map
#
#
# def query_acc(acc_ids):
#     id_pdb_map = _query(acc_ids, "ACC+ID")
#     return id_pdb_map
#
#
# def _query(input_names, tag):
#     assert isinstance(input_names, list)
#     query_line = " ".join(input_names)
#     url = 'https://www.uniprot.org/uploadlists/'
#
#     params = {'from': tag, 'to': 'PDB_ID', 'format': 'tab',
#               'query': query_line}
#     data = urllib.parse.urlencode(params)
#
#     data = data.encode('utf-8')
#     req = urllib.request.Request(url, data)
#     with urllib.request.urlopen(req) as f:
#         response = f.read()
#     parsed_response = response.decode('utf-8')
#     orig_id_pdb = parsed_response.strip().split("\n")
#     if len(orig_id_pdb) == 1: # Only header line
#         return []
#     # ['MSHA_SALTO\tmshA', 'NADE_STRCO\tnadE']
#     orig_id_pdb = orig_id_pdb[1:]
#     id_pdb_map = dict()
#     for line in orig_id_pdb:
#         input_id, pdb = line.split("\t")
#         id_pdb_map[input_id] = pdb
#     return id_pdb_map

# def convert_pdb_to_acc(pdb_list):
#     pdb_acc_map = uniprot_id_converter.convert("PDB_ID", "ACC", pdb_list)


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

