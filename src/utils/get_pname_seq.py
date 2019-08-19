"""
Obtain a pdb_code-seq map from sequence .fasta files
"""
import logging

from utils import generic, id_to_pdb

def parse_raw(seq_path):
    """
    Only uniprot
    """
    generic.quit_if_missing(seq_path)
    id_seq_map = dict()
    with open(seq_path, 'r') as file:
        current_seq = []
        for line in file:
            if line.startswith(">"):
                if current_seq:
                    seq = "".join(current_seq)
                    id_seq_map[desired_id] = seq
                    current_seq = []
                desired_id = line.strip().split("|")[1]
            else:
                current_seq.append(line.strip())
    if current_seq:
        id_seq_map[desired_id] = seq
    return id_seq_map
    pass


def parse(seq_path):
    """
    Figures out the sequence filetype automatically. Only distinguishes
    between uniprot and uniref now.
    """
    generic.quit_if_missing(seq_path)
    is_uniref = False
    is_valid = False
    with open(seq_path) as file:
        for line in file:
            if not line.startswith(">"):
                continue
            is_valid = True
            if line[1:].startswith("UniRef"):
                is_uniref = True
                break
            break
    if not is_valid:
        logging.error(f"Input Seq-file in {seq_path} is invalid, no header "
                      f"lines with > found.")
        raise Exception
    if is_uniref:
        pdb_seq_map = parse_uniref(seq_path)
    else:
        pdb_seq_map = parse_uniprot(seq_path)
    return pdb_seq_map


def parse_uniprot(seq_path):
    pdb_seq_map = _parse(seq_path, "uniprot")
    return pdb_seq_map


def parse_uniref(seq_path):
    pdb_seq_map = _parse(seq_path, "uniref")
    return pdb_seq_map


def _parse(seq_path, source):
    """
    Disable "uniref" for now, until we get the query to server for
    GENENAME=>PDB_ID working.
    """
    assert source in ("uniref", "uniprot")
    assert source in ("uniprot")
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
        except KeyError:
            # original_id does not have a corresponding pdb_code
            continue
        pdb_code = pdb_code.lower()
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
                current_seq.append(line.strip())
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


