import os
import pickle
import unittest

import numpy as np
import pandas as pd

from biopython_adapted import bio_interface
import preprocess
from tests.src import paths_test
from utils import generic
from meme_suite import meme_interface


class TestConverge(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.seq_path = paths_test.CONV_TEST_SEQ
        self.output = os.path.join(self.debug_folder, "output_1.pkl")
        self.binding_sites = paths_test.IONCOM_BINDING_SITES
        self.seed_seqs = paths_test.INPUT_CONV_SEED_SEQS
        self.ref_coverge_bind = paths_test.REF_CONV_WITH_BIND
        self.ref_converge_seed = paths_test.REF_CONV_WITH_SEED
        self.success = False

    def test_run_converge_bind(self):
        preprocess.run_converge(self.seq_path,
                                self.output,
                                binding_sites=self.binding_sites,
                                num_p=7,
                                storage_path=self.debug_folder)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_coverge_bind, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertEqual(len(ref_output), len(act_output))
        for pdb_info, (ref_df1, ref_df2, ref_df3) in ref_output.items():
            act_df1, act_df2, act_df3 = act_output[pdb_info]
            if ref_df1 is None:
                self.assertTrue(act_df1 is None)
            else:
                pd.testing.assert_frame_equal(ref_df1, act_df1)
            if ref_df2 is None:
                self.assertTrue(act_df2 is None)
            else:
                pd.testing.assert_frame_equal(ref_df2, act_df2)
            if ref_df3 is None:
                self.assertTrue(act_df3 is None)
            else:
                pd.testing.assert_frame_equal(ref_df3, act_df3)
        self.success = True

    def test_run_converge_seed(self):
        preprocess.run_converge(self.seq_path, self.output,
                                seed_seqs=self.seed_seqs, num_p=7,
                                storage_path=self.debug_folder)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_coverge_seed, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertEqual(len(ref_output), len(act_output))
        for pdb_info, (ref_df1, ref_df2, ref_df3) in ref_output.items():
            act_df1, act_df2, act_df3 = act_output[pdb_info]
            if ref_df1 is None:
                self.assertTrue(act_df1 is None)
            else:
                pd.testing.assert_frame_equal(ref_df1, act_df1)
            if ref_df2 is None:
                self.assertTrue(act_df2 is None)
            else:
                pd.testing.assert_frame_equal(ref_df2, act_df2)
            if ref_df3 is None:
                self.assertTrue(act_df3 is None)
            else:
                pd.testing.assert_frame_equal(ref_df3, act_df3)
        self.success = True

    def tearDown(self):
        if self.success:
            pass


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
            pass

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


class TestRunAll(unittest.TestCase):
    def setUp(self):
        self.debug_folder = generic.setup_debug_folder(paths_test.DEBUG)
        self.prosite_input = paths_test.PROSITE_EXTRACT
        self.ioncom_input = paths_test.IONCOM_EXTRACT
        self.output = os.path.join(self.debug_folder, "output_1.pkl")
        self.ref_meme_txt = paths_test.REF_MEME_TXT
        self.ref_prosite_meme = paths_test.REF_RUN_PROSITE_MEME
        self.ref_prosite_mast = paths_test.REF_RUN_PROSITE_MAST
        self.ref_ioncom_mast = paths_test.REF_RUN_IONCOM_MAST
        self.success = False

    def test_run_prosite_meme(self):
        preprocess.run_prosite_meme(extract_path=self.prosite_input,
                                    motif_len=13,
                                    output=self.output,
                                    num_p=7,
                                    storage_path=self.debug_folder)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_prosite_meme, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertEqual(len(ref_output), len(act_output))
        for pdb_info, (ref_df1, ref_df2, ref_df3) in ref_output.items():
            act_df1, act_df2, act_df3 = act_output[pdb_info]
            if ref_df1 is None:
                self.assertTrue(act_df1 is None)
            else:
                pd.testing.assert_frame_equal(ref_df1, act_df1)
            if ref_df2 is None:
                self.assertTrue(act_df2 is None)
            else:
                pd.testing.assert_frame_equal(ref_df2, act_df2)
            if ref_df3 is None:
                self.assertTrue(act_df3 is None)
            else:
                pd.testing.assert_frame_equal(ref_df3, act_df3)
        self.success = True

    def test_run_prosite_mast(self):
        preprocess.run_prosite_mast(extract_path=self.prosite_input,
                                    motif_len=13,
                                    ref_meme_txt=self.ref_meme_txt,
                                    output=self.output,
                                    storage_path=self.debug_folder)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_prosite_mast, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertEqual(len(ref_output), len(act_output))
        for pdb_info, (ref_df1, ref_df2, ref_df3) in ref_output.items():
            act_df1, act_df2, act_df3 = act_output[pdb_info]
            if ref_df1 is None:
                self.assertTrue(act_df1 is None)
            else:
                pd.testing.assert_frame_equal(ref_df1, act_df1)
            if ref_df2 is None:
                self.assertTrue(act_df2 is None)
            else:
                pd.testing.assert_frame_equal(ref_df2, act_df2)
            if ref_df3 is None:
                self.assertTrue(act_df3 is None)
            else:
                pd.testing.assert_frame_equal(ref_df3, act_df3)
        self.success = True

    def test_run_ioncom_mast(self):
        preprocess.run_ioncom_mast(extract_path=self.ioncom_input,
                                   motif_len=13,
                                   ref_meme_txt=self.ref_meme_txt,
                                   output=self.output,
                                   storage_path=self.debug_folder)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.ref_ioncom_mast, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertEqual(len(ref_output), len(act_output))
        for pdb_info, (ref_df1, ref_df2, ref_df3) in ref_output.items():
            act_df1, act_df2, act_df3 = act_output[pdb_info]
            if ref_df1 is None:
                self.assertTrue(act_df1 is None)
            else:
                pd.testing.assert_frame_equal(ref_df1, act_df1)
            if ref_df2 is None:
                self.assertTrue(act_df2 is None)
            else:
                pd.testing.assert_frame_equal(ref_df2, act_df2)
            if ref_df3 is None:
                self.assertTrue(act_df3 is None)
            else:
                pd.testing.assert_frame_equal(ref_df3, act_df3)
        self.success = True

    def tearDown(self):
        if self.success:
            pass


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
        preprocess.parse_extract_prosite(input_file=self.prosite_input,
                                         pname_cid_path=self.output)
        with open(self.output, 'rb') as file:
            act_output = pickle.load(file)
        with open(self.prosite_ref, 'rb') as file:
            ref_output = pickle.load(file)
        self.assertDictEqual(act_output, ref_output)
        self.success = True

    def test_ioncom(self):
        preprocess.parse_extract_ioncom(input_file=self.ioncom_input,
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

    def test_find_motif_meme(self):
        preprocess.find_motifs_meme(pname_cid_path=self.input_1,
                                    seq_file=self.seq_1,
                                    motif_len=13,
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

    def test_find_motif_mast(self):
        preprocess.find_motifs_mast(pname_cid_path=self.input_2,
                                    seq_file=self.seq_2,
                                    ref_meme_txt=self.meme_txt,
                                    motif_len=13,
                                    output=self.output,
                                    meme_folder=self.meme_folder)

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

