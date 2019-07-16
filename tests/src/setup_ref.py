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

Run setup_parse_extracts() before anything else that requires
pname_cid_map.
"""
import os

from preprocessing import preprocess
from tests.src import paths_test


# Preprocessing
def setup_preprocessing():
    setup_parse_extracts()
    setup_create_seq()
    setup_find_motif()
    setup_run_all()

def setup_parse_extracts():
    prosite_input = paths_test.PROSITE_EXTRACT
    prosite_ref = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
    preprocess.parse_extracts('prosite', prosite_input, prosite_ref)
    assert os.path.isfile(prosite_ref)

    ioncom_input = paths_test.IONCOM_EXTRACT
    ioncom_ref = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
    preprocess.parse_extracts('ioncom', ioncom_input, ioncom_ref)
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
    preprocess.find_motifs('meme',
                           pname_cid_path=input_1,
                           ref_meme_txt=None,
                           seq_file=seq_1,
                           output=output_1,
                           num_p=7)
    preprocess.find_motifs('mast',
                           pname_cid_path=input_2,
                           ref_meme_txt=paths_test.REF_MEME_TXT,
                           seq_file=seq_2,
                           output=output_2,
                           num_p=7)
    assert os.path.isfile(output_1)
    assert os.path.isfile(output_2)


def setup_run_all():
    input_1 = paths_test.PROSITE_EXTRACT
    input_2 = paths_test.IONCOM_EXTRACT
    output_1 = paths_test.REF_RUN_ALL_1
    output_2 = paths_test.REF_RUN_ALL_2
    output_3 = paths_test.REF_RUN_ALL_3
    ref_meme_txt = paths_test.REF_MEME_TXT
    preprocess.run_all(process='meme',
                       source='prosite',
                       num_p=7,
                       extract_path=input_1,
                       output=output_1,
                       delete_intermediate=True)
    preprocess.run_all(process='mast',
                       source='prosite',
                       extract_path=input_1,
                       ref_meme_txt=ref_meme_txt,
                       output=output_2,
                       delete_intermediate=True)
    preprocess.run_all(process='mast',
                       source='ioncom',
                       extract_path=input_2,
                       ref_meme_txt=ref_meme_txt,
                       output=output_3,
                       delete_intermediate=True)
    assert os.path.isfile(output_1)
    assert os.path.isfile(output_2)
    assert os.path.isfile(output_3)

def setup_meme_to_conv():
    input_meme_file = paths_test.ORIG_MEME
    ref_compositiom = paths_test.REF_CONV_COMP
    ref_output = paths_test.REF_CONV_OUTPUT
    preprocess.

if __name__ == "__main__":
    setup_preprocessing()