"""
Produces the reference output files for preprocess.py.

We have as input some form of information detailing the pname => cid mapping
and we wish to obtain a pname+cid => motif_pos at the end. First, we parse
the input information to obtain pname_cid_map. Together with a pre-downloaded
list of pdb files in pdb_folder, and a pre-computed meme.txt obtained from
running meme on the full ef-hand dataset or obtained elsewhere, this is fed
into parse_extracts() to obtain the desired output. This is then saved
separately for each source.

For prosite, the list of pdb codes representing the sequences are in
https://prosite.expasy.org/PS00018, copying the >>more list in PDB tab,
in the table. This forms pdb_list. We then take the html by
view_page_source in
https://prosite.expasy.org/cgi-bin/pdb/pdb_structure_list.cgi?src=PS00018,
and copy it into prosite_extract.txt.

For ioncom, we download the dataset from:
https://zhanglab.ccmb.med.umich.edu/IonCom/. This already provides us with
the pname+cid directly, so there's no need for a separate pdb_list.
"""
import pickle
import os
import shutil

import preprocess
from converge import conv_interface
from tests.src import paths_test
from meme_suite import meme_interface
from utils import generic, get_pname_seq, download_fasta_given_pdb
from biopython_adapted import bio_interface
from config import paths

def setup_all():
    # setup_parse_extracts()
    # setup_create_seq()
    # setup_find_motif()
    # setup_run_all()
    # setup_meme_suite_mast()
    # setup_meme_suite_meme()
    # setup_meme_to_conv()
    # setup_conv_to_meme()
    # setup_run_converge()
    # setup_get_pname_seq()
    pass

# todo: the one that broke is id_to_pdb, query for gene_name gave some weird
#  results apparently.

# todo: we can skip uniref for now, and rely entirely on uniprot. But note
#  that we also need to download pdb. So, test code need to change to calling
#  from the original list, instead of the new one.

def setup_get_pname_seq(to_run_uniref=False):
    uniprot_template_seq = paths_test.UNIPROT_SEQ
    if not os.path.isfile(uniprot_template_seq):
        #     This is how the uniprot_template_seq is generated. We want a seq
        #     for which the relevant .pdb files exist (so we don't have to
        #     download those just for testing), and the easiest way is to use
        #     the
        #     prosite_extract ones. But, this takes ~5min to run.
        prosite_input = paths.PROSITE_EXTRACT
        prosite_output_full = paths_test.TMP_FILE_TEMPLATE.format(1)
        preprocess.parse_extract_prosite(prosite_input, prosite_output_full)
        with open(prosite_output_full, 'rb') as file:
            pname_cid_map = pickle.load(file)
        pnames = list(pname_cid_map.keys())
        download_fasta_given_pdb.download(pnames, uniprot_template_seq)
    ref_uniprot_pname_seq = paths_test.REF_UNIPROT_PNAME_SEQ
    pname_seq_map = get_pname_seq.parse(uniprot_template_seq)
    with open(ref_uniprot_pname_seq, 'wb') as file:
        pickle.dump(pname_seq_map, file, -1)
    assert os.path.isfile(ref_uniprot_pname_seq)
    if to_run_uniref:
        # Uniref disabled for now, until we can get the parse() working for
        # uniref.
        uniref_seq = paths_test.UNIREF_SEQ
        ref_uniref_pname_seq = paths_test.REF_UNIREF_PNAME_SEQ
        pname_seq_map = get_pname_seq.parse(uniref_seq)
        with open(ref_uniref_pname_seq, 'wb') as file:
            pickle.dump(pname_seq_map, file, -1)
        assert os.path.isfile(ref_uniref_pname_seq)


# def setup_run_converge():
#     debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
#     input_seqs = paths_test.UNIPROT_SEQ
#     binding_sites = paths_test.IONCOM_BINDING_SITES
#     seed_seqs = paths_test.INPUT_CONV_SEED_SEQS
#     ref_conv_bind = paths_test.REF_CONV_WITH_BIND
#     ref_conv_seed = paths_test.REF_CONV_WITH_SEED
#
#     preprocess.run_converge(input_seqs,
#                             ref_conv_bind,
#                             binding_sites=binding_sites,
#                             num_p=7,
#                             storage_path=debug_folder)
    # preprocess.run_converge(input_seqs,
    #                         ref_conv_seed,
    #                         seed_seqs=seed_seqs,
    #                         num_p=7,
    #                         storage_path=debug_folder)
#     with open(ref_conv_seed, 'rb') as file:
#         print(pickle.load(file))
#     generic.quit_if_missing(ref_conv_bind)
#     generic.quit_if_missing(ref_conv_seed)


def setup_meme_suite_mast():
    debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
    mast_output_folder = os.path.join(debug_folder, "output_mast")
    input_seqfile = paths_test.MEME_TEST_SEQ
    input_memefile = paths_test.REF_MEME_TXT
    diagrams_output = paths_test.REF_MAST_DIAGRAMS

    meme_interface.run_mast(input_memefile, input_seqfile, mast_output_folder)
    mast_txt_path = os.path.join(mast_output_folder, "mast.txt")
    diagrams = bio_interface.parse_mast_file(mast_txt_path)
    with open(diagrams_output, 'wb') as file:
        pickle.dump(diagrams, file, -1)

def setup_meme_suite_meme():
    debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
    meme_output_folder = os.path.join(debug_folder, "output_meme")
    input_seqfile = paths_test.MEME_TEST_SEQ
    counts_output = paths_test.REF_MEME_COUNTS

    meme_interface.run_meme(input_seqfile, 30, meme_output_folder, num_p=7)
    meme_txt_path = os.path.join(meme_output_folder, "meme.txt")
    counts = bio_interface.parse_meme_file(meme_txt_path)
    with open(counts_output, 'wb') as file:
        pickle.dump(counts, file, -1)


def setup_parse_extracts():
    prosite_input = paths_test.PROSITE_EXTRACT
    prosite_ref = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
    preprocess.parse_extract_prosite(prosite_input, prosite_ref)
    assert os.path.isfile(prosite_ref)

    ioncom_input = paths_test.IONCOM_EXTRACT
    ioncom_ref = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
    preprocess.parse_extract_ioncom(ioncom_input, ioncom_ref)
    assert os.path.isfile(ioncom_ref)


def setup_create_seq():
    input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
    input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
    seq_1 = paths_test.REF_CREATE_SEQ_1
    seq_2 = paths_test.REF_CREATE_SEQ_2
    preprocess.create_seq(pname_cid_path=input_1,
                          pdb_folder=paths_test.PDB_FOLDER,
                          seq_path=seq_1)
    preprocess.create_seq(pname_cid_path=input_2,
                          pdb_folder=paths_test.PDB_FOLDER,
                          seq_path=seq_2)
    assert os.path.isfile(seq_1)
    assert os.path.isfile(seq_2)


def setup_find_motif():
    input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
    input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
    seq_1 = paths_test.REF_CREATE_SEQ_1
    seq_2 = paths_test.REF_CREATE_SEQ_2
    output_1 = paths_test.REF_FIND_MOTIF_1
    output_2 = paths_test.REF_FIND_MOTIF_2
    meme_folder = os.path.join(paths_test.DEBUG, 'meme_folder',
                               generic.get_timestamp())
    preprocess.find_motifs_meme(pname_cid_path=input_1,
                                seq_file=seq_1,
                                motif_len=13,
                                output=output_1,
                                meme_folder=meme_folder,
                                num_p=7)
    shutil.move(meme_folder, paths.TRASH)

    meme_folder = os.path.join(paths_test.DEBUG, 'meme_folder',
                               generic.get_timestamp())

    preprocess.find_motifs_mast(pname_cid_path=input_2,
                                seq_file=seq_2,
                                ref_meme_txt=paths_test.REF_MEME_TXT,
                                motif_len=13,
                                output=output_2,
                                meme_folder=meme_folder)
    shutil.move(meme_folder, paths.TRASH)
    assert os.path.isfile(output_1)
    assert os.path.isfile(output_2)


def setup_run_all():
    input_1 = paths_test.PROSITE_EXTRACT
    input_2 = paths_test.IONCOM_EXTRACT
    output_1 = paths_test.REF_RUN_PROSITE_MEME
    output_2 = paths_test.REF_RUN_PROSITE_MAST
    output_3 = paths_test.REF_RUN_IONCOM_MAST
    ref_meme_txt = paths_test.REF_MEME_TXT

    preprocess.run_prosite_meme(extract_path=input_1,
                                motif_len=13,
                                output=output_1,
                                num_p=7)
    preprocess.run_prosite_mast(extract_path=input_1,
                                motif_len=13,
                                ref_meme_txt=ref_meme_txt,
                                output=output_2)
    preprocess.run_ioncom_mast(extract_path=input_2,
                               motif_len=13,
                               ref_meme_txt=ref_meme_txt,
                               output=output_3)
    assert os.path.isfile(output_1)
    assert os.path.isfile(output_2)
    assert os.path.isfile(output_3)

def setup_meme_to_conv():
    meme_file_full = paths_test.ORIG_MEME_FOR_CONVERT
    ref_composition = paths_test.REF_CONV_COMPOSITION
    ref_matrix = paths_test.REF_CONV_MATRIX
    conv_interface.convert_meme_to_conv(meme_file_full, ref_composition,
                                        ref_matrix, minimal=False)
    assert os.path.isfile(ref_composition)
    assert os.path.isfile(ref_matrix)


def setup_meme_to_conv_minimal():
    """
    Technically not required at this point, since earlier 2 setup funcs has
    already made the necessary files. It also overwrites existing ref using a
    different method.
    """
    meme_file = paths_test.REF_MEME_FROM_CONV
    ref_composition = paths_test.REF_CONV_COMPOSITION
    ref_matrix = paths_test.REF_CONV_MATRIX
    conv_interface.convert_meme_to_conv(meme_file, ref_composition,
                                        ref_matrix, minimal=True)
    assert os.path.isfile(ref_composition)
    assert os.path.isfile(ref_matrix)

def setup_conv_to_meme():
    input_composition = paths_test.REF_CONV_COMPOSITION
    input_matrix = paths_test.REF_CONV_MATRIX
    ref_meme = paths_test.REF_MEME_FROM_CONV
    conv_interface.convert_conv_to_meme(input_matrix, input_composition,
                                        ref_meme)


if __name__ == "__main__":
    setup_all()