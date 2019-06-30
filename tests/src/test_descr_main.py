import logging
import unittest
import os
import pickle
import shutil

from descr import descr_main, loaders
from utils import extract_parser, motif_finder
import config
import pandas as pd

# todo: loaders.load_pdb_info also need to test


class TestDescrCalculate(unittest.TestCase):
    def setUp(self):
        ref_ = os.path.join(config.ROOT, 'tests', 'data', 'ref')

        self.pdb_folder_path = os.path.join(config.ROOT, 'data', 'input',
                                            'pdb_files')

        p_motif_pos_dir = os.path.join(ref_, 'motif_pos_prosite.pkl')
        i_motif_pos = os.path.join(ref_, 'motif_pos_ioncom.pkl')
        p_mast_motif_pos = os.path.join(ref_, 'motif_pos_meme_prosite.pkl')

        with open(p_motif_pos_dir, 'rb') as file:
            self.input_prosite_meme = pickle.load(file)

        with open(p_mast_motif_pos, 'rb') as file:
            self.input_prosite_mast = pickle.load(file)

        with open(i_motif_pos, 'rb') as file:
            self.input_ioncom_mast = pickle.load(file)

        ref_prosite_meme_dir = os.path.join(ref_, 'descr_prosite_meme.pkl')
        ref_prosite_mast_dir = os.path.join(ref_, 'descr_prosite_mast.pkl')
        ref_ioncom_mast_dir = os.path.join(ref_, 'descr_ioncom_mast.pkl')

        with open(ref_prosite_meme_dir, 'rb') as file:
            self.ref_prosite_meme = pickle.load(file)
        with open(ref_prosite_mast_dir, 'rb') as file:
            self.ref_prosite_mast = pickle.load(file)
        with open(ref_ioncom_mast_dir, 'rb') as file:
            self.ref_ioncom_mast = pickle.load(file)

    def test_prosite_meme(self):
        pdb_info = loaders.load_pdb_info(self.input_prosite_meme,
                                         pdb_dir=self.pdb_folder_path)
        descrs = descr_main.calculate(pdb_info)
        pd.testing.assert_frame_equal(descrs, self.ref_prosite_meme)

    def test_prosite_mast(self):
        pdb_info = loaders.load_pdb_info(self.input_prosite_mast,
                                                  pdb_dir=self.pdb_folder_path)
        descrs = descr_main.calculate(pdb_info)
        pd.testing.assert_frame_equal(descrs, self.ref_prosite_mast)

    def test_ioncom_mast(self):
        pdb_info = loaders.load_pdb_info(self.input_ioncom_mast,
                                                  pdb_dir=self.pdb_folder_path)
        descrs = descr_main.calculate(pdb_info)
        pd.testing.assert_frame_equal(descrs, self.ref_ioncom_mast)
