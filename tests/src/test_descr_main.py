# todo: need to test loaders.load_pdb_info

# import logging
import unittest
import os
import pickle

import pandas as pd

from descr import descr_main, loaders
import config

class TestDescrCalculate(unittest.TestCase):
    def setUp(self):
        ref_ = os.path.join(config.ROOT, 'tests', 'data', 'ref')

        self.pdb_folder_path = os.path.join(config.ROOT, 'data', 'input',
                                            'pdb_files')

        self.p_motif_pos_dir = os.path.join(ref_, 'motif_pos_prosite.pkl')
        self.i_motif_pos = os.path.join(ref_, 'motif_pos_ioncom.pkl')
        self.p_mast_motif_pos = os.path.join(ref_, 'motif_pos_meme_prosite.pkl')

        self.ref_prosite_meme_dir = os.path.join(ref_, 'descr_prosite_meme.pkl')
        self.ref_prosite_mast_dir = os.path.join(ref_, 'descr_prosite_mast.pkl')
        self.ref_ioncom_mast_dir = os.path.join(ref_, 'descr_ioncom_mast.pkl')

    def test_prosite_meme(self):
        with open(self.p_motif_pos_dir, 'rb') as file:
            input_prosite_meme = pickle.load(file)
        with open(self.ref_prosite_meme_dir, 'rb') as file:
            ref_prosite_meme = pickle.load(file)

        pdb_info = loaders.load_pdb_info(input_prosite_meme,
                                         pdb_dir=self.pdb_folder_path)
        descrs = descr_main.calculate(pdb_info)
        pd.testing.assert_frame_equal(descrs, ref_prosite_meme)

    def test_prosite_mast(self):
        with open(self.p_mast_motif_pos, 'rb') as file:
            input_prosite_mast = pickle.load(file)
        with open(self.ref_prosite_mast_dir, 'rb') as file:
            ref_prosite_mast = pickle.load(file)

        pdb_info = loaders.load_pdb_info(input_prosite_mast,
                                         pdb_dir=self.pdb_folder_path)
        descrs = descr_main.calculate(pdb_info)
        pd.testing.assert_frame_equal(descrs, ref_prosite_mast)

    def test_ioncom_mast(self):
        with open(self.i_motif_pos, 'rb') as file:
            input_ioncom_mast = pickle.load(file)
        with open(self.ref_ioncom_mast_dir, 'rb') as file:
            ref_ioncom_mast = pickle.load(file)

        pdb_info = loaders.load_pdb_info(input_ioncom_mast,
                                         pdb_dir=self.pdb_folder_path)
        descrs = descr_main.calculate(pdb_info)
        pd.testing.assert_frame_equal(descrs, ref_ioncom_mast)
