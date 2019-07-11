import logging
import os

# pylint: disable=invalid-name
def extract_sequence(pdb_folder, pname_cid_map, AA3_to_AA1):
    """
    Extracts the seq for each pname-cid from their .pdb file.
    """
    seqs = []
    for pname in pname_cid_map.keys():
        cid = pname_cid_map[pname]
        assert len(cid) == 1
        seq = [f">{pname}\n"]
        pdb_filepath = os.path.join(pdb_folder, pname + ".pdb")
        with open(pdb_filepath, 'r') as pdb_file:
            current_atom_count = 1
            for line in pdb_file:
                if line.startswith("MODEL        2"):
                    break
                if not line.startswith("ATOM"):
                    continue
                if line[21] != cid:
                    continue
                if line[26] != " ":  # altloc
                    continue
                sno = int(line[22:26].strip())
                if sno < current_atom_count:
                    continue
                elif sno > current_atom_count:
                    num_spacer_res = sno - current_atom_count
                    seq.append("X" * num_spacer_res)
                    current_atom_count = sno
                res = line[17:20].strip()
                AA_single_letter = AA3_to_AA1[res]
                seq.append(f"{AA_single_letter}")
                current_atom_count += 1
            seq.append("\n")
        seq = "".join(seq)
        seqs.append(seq)
    return seqs
# pylint: enable=invalid-name

def test_fasta_match_pdb(fasta_fname, pdb_folder, seq_cid_map, AA3_to_AA1):
    pname = None
    with open(fasta_fname, 'r') as f_file:
        for f_line in f_file:
            if f_line.startswith(">"):
                pname = f_line[1:5]
                continue
            if pname:
                pdb_path = os.path.join(pdb_folder, pname + '.pdb')
                cid = seq_cid_map[pname]
                with open(pdb_path, 'r') as p_file:
                    for p_line in p_file:
                        if not p_line.startswith("ATOM"):
                            continue
                        if p_line[21] != cid:
                            continue
                        if p_line[26] != " ":  # altloc
                            continue
                        sno = int(p_line[22:26].strip())
                        res = p_line[17:20]
                        if sno < 1:
                            continue
                        if AA3_to_AA1[res] != f_line[sno - 1]:
                            logging.error(
                                f"Mismatch between seq in pdb file and in "
                                f"extracted fasta file. pdb_line: <{p_line}>, "
                                f"fasta_line: <{f_line}>")
                            raise Exception