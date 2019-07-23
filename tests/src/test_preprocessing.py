import os
import pickle
import unittest
import shutil

import numpy as np

import preprocess
from tests.src import paths_test
from utils import generic
from meme_suite import meme_interface
from biopython_adapted import bio_interface


class TestMemeSuiteMast(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.input_seqfile = paths_test.MEME_TEST_SEQ
        self.input_memefile = paths_test.REF_MEME_TXT
        self.output = os.path.join(self.debug_folder, "output_meme")
        self.ref_output = paths_test.REF_MAST_DIAGRAMS
        self.success = False

    def test_mast(self):
        meme_interface.run_mast(self.input_memefile, self.input_seqfile,
                                self.output)
        mast_txt_path = os.path.join(self.output, "mast.txt")
        act_diagrams = bio_interface.parse_mast_file(mast_txt_path)
        with open(self.ref_output, 'rb') as file:
            ref_diagrams = pickle.load(file)
        self.assertDictEqual(act_diagrams, ref_diagrams)
        self.success = True

    def tearDown(self):
        if self.success:
            pass  # shutil.rmtree(self.debug_folder)

class TestMemeSuiteMeme(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.input_seqfile = paths_test.MEME_TEST_SEQ
        self.output = os.path.join(self.debug_folder, "output_meme")
        self.ref_output = paths_test.REF_MEME_COUNTS

        self.success = False

    def test_meme(self):
        meme_interface.run_meme(self.input_seqfile, 30, self.output, num_p=7)
        meme_txt_path = os.path.join(self.output, "meme.txt")
        act_counts = bio_interface.parse_meme_file(meme_txt_path)
        with open(self.ref_output, 'rb') as file:
            ref_counts = pickle.load(file)
        np.testing.assert_array_equal(act_counts, ref_counts)
        self.success = True

    def tearDown(self):
        if self.success:
            pass
            # shutil.rmtree(self.debug_folder)


class TestFindCid(unittest.TestCase):
    pass


class TestRunAll(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.prosite_input = paths_test.PROSITE_EXTRACT
        self.ioncom_input = paths_test.IONCOM_EXTRACT
        self.output = os.path.join(self.debug_folder, "output_1.pkl")
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
                     output=self.output, storage_path=self.debug_folder)
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
                     output=self.output, storage_path=self.debug_folder)
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
                     output=self.output, storage_path=self.debug_folder)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_output_3, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def tearDown(self):
        if self.success:
            pass
            # shutil.rmtree(self.debug_folder)


class TestParseExtract(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.prosite_input = paths_test.PROSITE_EXTRACT
        self.ioncom_input = paths_test.IONCOM_EXTRACT
        self.output = os.path.join(self.debug_folder, "output_1.pkl")

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
            pass
            # shutil.rmtree(self.debug_folder)


class TestCreateSeq(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
        self.input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
        self.pdb_folder = paths_test.PDB_FOLDER
        self.output = os.path.join(self.debug_folder, "output_1.pkl")

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
            pass
            # shutil.rmtree(self.debug_folder)


class TestFindMotif(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.input_1 = paths_test.REF_PROSITE_EXTRACT_PNAME_CID
        self.input_2 = paths_test.REF_IONCOM_EXTRACT_PNAME_CID
        self.seq_1 = paths_test.REF_CREATE_SEQ_1
        self.seq_2 = paths_test.REF_CREATE_SEQ_2
        self.pdb_folder = paths_test.PDB_FOLDER
        self.output = os.path.join(self.debug_folder, "output_1.pkl")
        self.meme_folder = os.path.join(self.debug_folder, "meme_folder")
        self.meme_txt = paths_test.REF_MEME_TXT

        self.ref_output_1 = paths_test.REF_FIND_MOTIF_1
        self.ref_output_2 = paths_test.REF_FIND_MOTIF_2

        self.success = False

    def test_find_motif_1(self):
        preprocess.find_motifs('meme', 13,
                               pname_cid_path=self.input_1,
                               ref_meme_txt=None,
                               seq_file=self.seq_1,
                               output=self.output,
                               meme_folder=self.meme_folder,
                               num_p=7)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_output_1, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertIsInstance(act_output, dict)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def test_find_motif_2(self):
        preprocess.find_motifs('mast', 13,
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
            pass
            # shutil.rmtree(self.debug_folder)

