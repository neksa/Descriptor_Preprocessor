"""
Obtain a pname/acc-seq map from sequence .fasta files
"""
from utils import generic, id_to_pdb

def parse_uniprot(seq_path):
    pdb_seq_map = _parse(seq_path, "uniprot")
    return pdb_seq_map


def parse_uniref(seq_path):
    pdb_seq_map = _parse(seq_path, "uniref")
    return pdb_seq_map


def _parse(seq_path, source):
    assert source in ("uniref", "uniprot")
    if source == "uniref":
        extractor = _extract_genename
        query_fn = id_to_pdb.query_genename
    else:
        extractor = _extract_accid
        query_fn = id_to_pdb.query_acc
    id_seq_map = _get_id_seq_map(seq_path, extractor)
    orig_ids = list(id_seq_map.keys())
    id_pdb_map = query_fn(orig_ids)
    pdb_seq_map = dict()
    for orig_id, seq in id_seq_map.items():
        try:
            pdb_code = id_pdb_map[orig_id]
        except IndexError:
            # original_id does not have a corresponding pdb_code
            continue
        pdb_seq_map[pdb_code] = seq
    return pdb_seq_map


def _get_id_seq_map(seq_path, line_extractor):
    generic.quit_if_missing(seq_path)
    id_seq_map = dict()
    with open(seq_path, 'r') as file:
        current_seq = []
        for line in file:
            if line.startswith(">"):
                desired_id = line_extractor(line)
                if current_seq:
                    seq = "".join(current_seq)
                    id_seq_map[desired_id] = seq
                    current_seq = []
            else:
                current_seq.append(line)
    if current_seq:
        id_seq_map[desired_id] = seq
    return id_seq_map


def _extract_genename(raw_line):
    split_line = raw_line.strip().split("RepID=")
    assert len(split_line) == 2
    gene_name = split_line[1]
    return gene_name


def _extract_accid(raw_line):
    split_line = raw_line.strip().split("|")
    acc_id = split_line[1]
    return acc_id


