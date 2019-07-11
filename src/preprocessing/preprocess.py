import pickle
import os
import logging

from config import paths
from preprocessing import extract_parser, download_pdb_files, \
    create_seq_file, filter_seqs, make_conv_seed_seqs, motif_finder, \
    prosite_pdb_list
from utils import generic

def run_all(process='meme', source='prosite', num_p=7, extract_path=None,
            pname_cid_path=paths.PNAME_CID, pdb_folder=paths.PDB_FOLDER,
            seq_path=paths.FULL_SEQS, ref_meme_txt=paths.REF_MEME_TXT,
            mast_meme_folder=paths.MEME_MAST_FOLDER, output=paths.MOTIF_POS):
    assert isinstance(num_p, int)
    assert num_p >= 1
    assert process in ('meme', 'mast')
    assert source in ('prosite', 'ioncom')

    parse_extracts(source, extract_path, pname_cid_path)
    download_pdb(pname_cid_path, pdb_folder)
    trim_pnames_based_on_pdb(pname_cid_path, pdb_folder)
    create_seq(pname_cid_path, pdb_folder, seq_path)

    filter_seq_file(seq_path)
    # create_conv_seed_seqs()
    find_motifs(process, pname_cid_path, ref_meme_txt, mast_meme_folder,
                seq_path, output,
                num_p)


def parse_extracts(source='prosite', input_file_path=None,
                   pname_cid_path=paths.PNAME_CID):
    assert source in ('prosite', 'ioncom')
    if input_file_path == None:
        if source == 'prosite':
            input_file_path = paths.PROSITE_EXTRACT
        else:
            input_file_path = paths.IONCOM_EXTRACT
    generic.quit_if_missing(input_file_path)
    if source == 'prosite':
        pname_cid_map = extract_parser.parse_prosite(input_file_path,
                                                     prosite_pdb_list.pdb_list)
    else:
        pname_cid_map = extract_parser.parse_ioncom(input_file_path)
    _test_seq_cid_map(pname_cid_map)
    generic.warn_if_exist(pname_cid_path)
    with open(pname_cid_path, 'wb') as file:
        pickle.dump(pname_cid_map, file, -1)


def download_pdb(pname_cid_path=paths.PNAME_CID, pdb_folder=paths.PDB_FOLDER,
                 replace_existing=False):
    if replace_existing:
        logging.warning(f"Replacing .pdb files in PDB_folder <{pdb_folder}>. "
                        f"Downloading from scratch takes a while.")
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    if not os.path.isdir(pdb_folder):
        logging.warning(f"PDB_folder in <{pdb_folder}> not found. Downloading "
                        f"from scratch takes a while.")
        os.mkdir(pdb_folder)
    download_pdb_files.download(pname_cid_map, pdb_folder, False)


def trim_pnames_based_on_pdb(pname_cid_path=paths.PNAME_CID,
                             pdb_folder=paths.PDB_FOLDER):
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    download_pdb_files.trim_pname_cid(pname_cid_map, pdb_folder)
    _test_seq_cid_map(pname_cid_map)
    with open(pname_cid_path, 'wb') as file:
        pickle.dump(pname_cid_map, file, -1)
    generic.quit_if_missing(pname_cid_path)


def create_seq(pname_cid_path=paths.PNAME_CID,
               pdb_folder=paths.PDB_FOLDER, seq_path=paths.FULL_SEQS):
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    seqs = create_seq_file.extract_sequence(pdb_folder, pname_cid_map,
                                            AA3_to_AA1=generic.AA3_to_AA1)
    with open(seq_path, 'w') as file:
        file.writelines(seqs)
    create_seq_file.test_fasta_match_pdb(seq_path, pdb_folder,
                                         pname_cid_map, generic.AA3_to_AA1)


def filter_seq_file(seq_path=paths.FULL_SEQS):
    generic.quit_if_missing(seq_path)
    filter_seqs.delete_short_seqs(seq_path, threshold=30)


def find_motifs(process,
                pname_cid_path=paths.PNAME_CID,
                ref_meme_txt=paths.REF_MEME_TXT,
                mast_meme_folder=paths.MEME_MAST_FOLDER,
                seq_file=paths.FULL_SEQS,
                output=paths.MOTIF_POS,
                num_p=1):
    assert process in ('mast', 'meme')
    generic.quit_if_missing(pname_cid_path)
    with open(pname_cid_path, 'rb') as file:
        pname_cid_map = pickle.load(file)
    _test_seq_cid_map(pname_cid_map)
    motif_pos = motif_finder.find(pname_cid_map,
                                  process=process,
                                  num_p=num_p,
                                  ref_meme_txt=ref_meme_txt,
                                  mast_meme_folder=mast_meme_folder,
                                  seq_file=seq_file)
    generic.warn_if_exist(output)
    with open(output, 'wb') as file:
        pickle.dump(motif_pos, file, -1)


def create_conv_seed_seqs(binding_site_path=paths.IONCOM_BINDING_SITES,
                          seed_seq_path=paths.CONV_SEED_SEQS):
    generic.quit_if_missing(binding_site_path)
    generic.quit_if_missing(seed_seq_path)
    make_conv_seed_seqs.make(binding_site_path, seed_seq_path)


def _test_seq_cid_map(seq_cid_map):
    for pname, cid in seq_cid_map.items():
        assert isinstance(pname, str)
        assert isinstance(cid, str)
        assert len(pname) == 4
        assert len(cid) == 1
    return True



# def plot_labelled_motifs(seq_motif_map, fasta_filename):
#     """
#     This debugging function displays the sequence logo of the identified
#     motifs.
#     """
#     motifs = _get_labelled_motifs(seq_motif_map, fasta_filename)
#     plt.figure()
#     converted_motifs = []
#     extracted_motifs = np.array([list(i) for i in motifs])
#     for AA_per_pos in extracted_motifs.T:
#         AA, counts = np.unique(AA_per_pos, return_counts=True)
#         sorted_i = np.argsort(counts)
#         counts = counts[sorted_i]
#         AA = AA[sorted_i]
#         probs = counts / np.sum(counts)
#         motif_this_pos = []
#         for aa, prob in zip(AA, probs):
#             motif_this_pos.append([aa, prob])
#         converted_motifs.append(motif_this_pos)
#     seq_logo.Logo(converted_motifs, -1, convert_AA3=False)
#
#
# def _get_labelled_motifs(seq_motif_map, fasta_filename):
#     labelled_motifs = []
#     pname = None
#     with open(fasta_filename, 'r') as file:
#         for line in file:
#             if line.startswith('>'):
#                 pname = line[1:5]
#                 continue
#             if pname in seq_motif_map:
#                 motif_pos = seq_motif_map[pname]
#             else:
#                 #   Might have been dropped because of gapped motif, etc
#                 pname = None
#                 continue
#             for pos in motif_pos:
#                 if len(line) <= pos + 13:
#                     logging.error(
#                         f"Fasta seq is shorter than pos+13, for pos in "
#                         f"motif_pos. Fasta_seq: <{line}>, "
#                         f"motif_pos: <{motif_pos}>, illegal pos: <{pos}>.")
#                     raise Exception
#                 labelled_motifs.append(line[pos - 1:pos + 12])
#     return labelled_motifs