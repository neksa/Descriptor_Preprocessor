import re

def make(ioncom_dir, seed_seq_path):
    seq_binding_sites = ioncom_parser(ioncom_dir)
    profiles = split_into_profiles(seq_binding_sites)
    profiles = sorted(profiles)
    write_seed_seq_file(profiles, seed_seq_path)

def ioncom_parser(fname):
    # Input: Ioncom binding site file
    # Output: list of 30-len profiles
    with open(fname, 'r') as file:
        lines = file.readlines()
    # Split into seq + binding sites
    seq_binding_sites = []
    seq = None
    binding_site = None
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if seq is not None:
                assert len(seq) == len(binding_site)
                seq_binding_sites.append([seq, binding_site])
                seq = None
            else:
                continue
        elif re.match('[A-Z]', line):
            if len(line) > 30:
                seq = line
        else:
            binding_site = line
    if seq:
        seq_binding_sites.append([seq, binding_site])
    return seq_binding_sites

def split_into_profiles(seq_binding_sites):
    # Input: List of seq + binding sites
    # Output: Split profiles, as a set, so screened for duplicates
    profiles = set()
    for seq, binding_site in seq_binding_sites:
        length = len(seq)
        loc_of_matches = [i.span()[0] for i in re.finditer("1", binding_site)]
        for index in loc_of_matches:
            if index < 15:
                profile = seq[:30]
            elif index + 15 >= length:
                profile = seq[-30:]
            else:
                profile = seq[index-15:index+15]
            assert len(profile) == 30
            profiles.add(profile)
    return profiles

def write_seed_seq_file(profiles, output_fname):
    with open(output_fname, 'w') as file:
        file.write(">SEED_SEQ\n")
        for i, profile in enumerate(profiles):
            # if i == 30:
            #     break
            file.write(profile)
    return

from config import paths
ioncom_dir = paths.ROOT + '/data/input/ioncom/allsulfate.txt'
seed_seq_path = paths.ROOT + '/data/seed_seqs.fasta'
make_seed_seqs(ioncom_dir, seed_seq_path)