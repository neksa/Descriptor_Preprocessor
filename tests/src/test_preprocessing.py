import os
import pickle
import unittest

from preprocessing import preprocess
from tests.src import paths_test

class TestFindCid(unittest.TestCase):
    pass


class TestRunAll(unittest.TestCase):
    def setUp(self):
        self.prosite_input = paths_test.PROSITE_EXTRACT
        self.ioncom_input = paths_test.IONCOM_EXTRACT
        self.output = paths_test.TMP_FILE
        self.ref_meme_txt = paths_test.REF_MEME_TXT
        self.ref_output_1 = paths_test.REF_RUN_ALL_1
        self.ref_output_2 = paths_test.REF_RUN_ALL_2
        self.ref_output_3 = paths_test.REF_RUN_ALL_3

        self.success = False

    def test_run_all_1(self):
        preprocess.run_all(process='meme',
                           source='prosite',
                           num_p=7,
                           extract_path=self.prosite_input,
                           output=self.output,
                           delete_intermediate=True)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_output_1, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def test_run_all_2(self):
        preprocess.run_all(process='mast',
                           source='prosite',
                           extract_path=self.prosite_input,
                           ref_meme_txt=self.ref_meme_txt,
                           output=self.output,
                           delete_intermediate = True)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_output_2, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def test_run_all_3(self):
        preprocess.run_all(process='mast',
                           source='ioncom',
                           extract_path=self.ioncom_input,
                           ref_meme_txt=self.ref_meme_txt,
                           output=self.output,
                           delete_intermediate=True)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_output_3, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def tearDown(self):
        if self.success:
            os.remove(self.output)


class TestParseExtract(unittest.TestCase):
    def setUp(self):
        self.prosite_input = paths_test.PROSITE_EXTRACT
        self.ioncom_input = paths_test.IONCOM_EXTRACT
        self.output = paths_test.TMP_FILE

        self.prosite_ref = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
        self.ioncom_ref = paths_test.REF_IONCOM_EXTRACT_PNAME_CID

        self.success = False

    def test_prosite(self):
        preprocess.parse_extracts(source='prosite',
                                  input_file_path=self.prosite_input,
                                  pname_cid_path=self.output)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.prosite_ref, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def test_ioncom(self):
        preprocess.parse_extracts(source='ioncom',
                                  input_file_path=self.ioncom_input,
                                  pname_cid_path=self.output)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ioncom_ref, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def tearDown(self):
        if self.success:
            os.remove(self.output)


class TestCreateSeq(unittest.TestCase):
    def setUp(self):
        self.input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
        self.input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
        self.pdb_folder = paths_test.PDB_FOLDER
        self.output = paths_test.TMP_FILE

        self.ref_output_1 = paths_test.REF_CREATE_SEQ_1
        self.ref_output_2 = paths_test.REF_CREATE_SEQ_2

        self.success = False

    def test_create_seq_1(self):
        preprocess.create_seq(pname_cid_path=self.input_1,
                              pdb_folder=self.pdb_folder,
                              seq_path=self.output)
        with open(self.output, 'r') as file:
            act_lines = file.readlines()
        with open(self.ref_output_1, 'r') as file:
            ref_lines = file.readlines()
        self.assertListEqual(act_lines, ref_lines)
        self.success = True

    def test_create_seq_2(self):
        preprocess.create_seq(pname_cid_path=self.input_2,
                              pdb_folder=self.pdb_folder,
                              seq_path=self.output)
        with open(self.output, 'r') as file:
            act_lines = file.readlines()
        with open(self.ref_output_2, 'r') as file:
            ref_lines = file.readlines()
        self.assertListEqual(act_lines, ref_lines)
        self.success = True

    def tearDown(self):
        if self.success:
            os.remove(self.output)


class TestFindMotif(unittest.TestCase):
    def setUp(self):
        self.input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
        self.input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
        self.seq_1 = paths_test.REF_CREATE_SEQ_1
        self.seq_2 = paths_test.REF_CREATE_SEQ_2
        self.pdb_folder = paths_test.PDB_FOLDER
        self.output = paths_test.TMP_FILE
        self.meme_txt = paths_test.REF_MEME_TXT

        self.ref_output_1 = paths_test.REF_FIND_MOTIF_1
        self.ref_output_2 = paths_test.REF_FIND_MOTIF_2

        self.success = False

    def test_find_motif_1(self):
        preprocess.find_motifs('meme',
                               pname_cid_path=self.input_1,
                               ref_meme_txt=None,
                               seq_file=self.seq_1,
                               output=self.output,
                               num_p=7)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_output_1, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertIsInstance(act_output, dict)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def test_find_motif_2(self):
        preprocess.find_motifs('mast',
                               pname_cid_path=self.input_2,
                               ref_meme_txt=self.meme_txt,
                               seq_file=self.seq_2,
                               output=self.output,
                               num_p=7)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_output_2, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertIsInstance(act_output, dict)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def tearDown(self):
        if self.success:
            os.remove(self.output)

