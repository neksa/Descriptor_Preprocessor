import logging
import unittest
import os
import pickle
import shutil

from utils import extract_parser, motif_finder
import config

class TestIoncomExtract(unittest.TestCase):
    def setUp(self):
        ref_pname_cid_path = os.path.join(config.ROOT, 'tests', 'data',
                                          'ref', 'pname_cid_ioncom.pkl')
        with open(ref_pname_cid_path, 'rb') as file:
            self.ref_pname_cid = pickle.load(file)

    def test_parse_ioncom(self):
        ioncom_input_path = os.path.join(config.ROOT, 'tests', 'data',
                                         'input', 'ioncom.txt')
        self.assertTrue(os.path.isfile(ioncom_input_path))
        act_pname_cid_map = extract_parser.parse_ioncom(ioncom_input_path)
        self.assertIsInstance(act_pname_cid_map, dict)
        self.assertDictEqual(self.ref_pname_cid, act_pname_cid_map)

class TestPrositeExtract(unittest.TestCase):
    def setUp(self):
        ref_pname_cid_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                          'pname_cid_prosite.pkl')
        with open(ref_pname_cid_path, 'rb') as file:
            self.ref_pname_cid = pickle.load(file)

    def test_parse_prosite(self):
        prosite_input = os.path.join(config.ROOT, 'tests', 'data', 'input',
                                     'prosite_extract.txt')
        self.assertTrue(os.path.isfile(prosite_input))
        act_pname_cid_map = extract_parser.parse_prosite(
            prosite_input, config.prosite_pdb_list)
        self.assertIsInstance(act_pname_cid_map, dict)
        self.assertDictEqual(self.ref_pname_cid, act_pname_cid_map)

class TestMotifFinder(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = os.path.join(config.ROOT, 'tests', 'data',
                                    'tmp_MotifFinder')
        if os.path.isdir(self.tmp_dir):
            logging.warning(f"Test {self.tmp_dir} folder not deleted "
                            f"previously. Deleting.")
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)
        os.mkdir(os.path.join(self.tmp_dir, 'input'))

        self.store_dir = os.path.join(self.tmp_dir, 'store')
        os.mkdir(self.store_dir)

        ref_folder = os.path.join(config.ROOT, 'tests', 'data', 'ref')

        pname_cid_ioncom_path = os.path.join(
            ref_folder, 'pname_cid_ioncom.pkl')
        self.assertTrue(os.path.isfile(pname_cid_ioncom_path))
        with open(pname_cid_ioncom_path, 'rb') as file:
            self.pname_cid_ioncom = pickle.load(file)

        pname_cid_prosite_path = os.path.join(
            ref_folder, 'pname_cid_prosite.pkl')
        self.assertTrue(os.path.isfile(pname_cid_prosite_path))
        with open(pname_cid_prosite_path, 'rb') as file:
            self.pname_cid_prosite = pickle.load(file)

        self.pdb_folder = os.path.join(config.ROOT, 'data', 'input',
                                       'pdb_files')

        motif_pos_prosite_path = os.path.join(
            ref_folder, 'motif_pos_prosite.pkl')
        with open(motif_pos_prosite_path, 'rb') as file:
            self.ref_motif_pos_prosite = pickle.load(file)

        motif_pos_ioncom_path = os.path.join(
            ref_folder, 'motif_pos_ioncom.pkl')
        with open(motif_pos_ioncom_path, 'rb') as file:
            self.ref_motif_pos_ioncom = pickle.load(file)

        self.ref_meme_path = os.path.join(ref_folder, 'meme.txt')
        self.success_1 = False
        self.success_2 = False
        self.success_3 = False

    def test_find_motif_pos_prosite(self):
        act_motif_pos = motif_finder.find_motif_pos(
            self.pname_cid_prosite, self.pdb_folder, process='meme',
            store_dir=self.store_dir, replace_existing=True)
        self.assertIsInstance(act_motif_pos, dict)
        self.assertDictEqual(act_motif_pos, self.ref_motif_pos_prosite)
        self.success_1 = True

    def test_find_motif_pos_ioncom(self):
        act_motif_pos = motif_finder.find_motif_pos(
            self.pname_cid_ioncom, self.pdb_folder, process='mast',
            store_dir=self.store_dir, replace_existing=True,
            ref_meme_txt=self.ref_meme_path)
        self.assertIsInstance(act_motif_pos, dict)
        self.assertDictEqual(act_motif_pos, self.ref_motif_pos_ioncom)
        self.success_2 = True

    def test_find_motif_pos_prosite_no_pdb_available(self):
        empty_pdb_folder = os.path.join(self.tmp_dir, 'empty_pdb_files')
        if os.path.isdir(empty_pdb_folder):
            logging.error(f"Why is empty_pdb_folder in <{empty_pdb_folder}> "
                          f"not deleted?")
            shutil.rmtree(empty_pdb_folder)
        os.mkdir(empty_pdb_folder)
        shortened_pname_cid = dict()
        for i, key in enumerate(self.pname_cid_prosite.keys()):
            if i == 4:
                break
            shortened_pname_cid[key] = self.pname_cid_prosite[key]
        motif_finder.find_motif_pos(shortened_pname_cid,
                                    empty_pdb_folder,
                                    process='meme',
                                    store_dir=self.store_dir,
                                    replace_existing=True,
                                    delete_intermediate_store=False)
        self.assertEqual(len(list(os.listdir(empty_pdb_folder))), 4)
        for file in os.listdir(empty_pdb_folder):
            path = os.path.join(empty_pdb_folder, file)
            self.assertTrue(os.path.isfile(path))
        self.success_3 = True

    def tearDown(self):
        if not os.path.isdir(self.tmp_dir):
            logging.warning(
                f"test_store_dir <{self.tmp_dir}> does not exists when "
                f"running TearDown. Check to see why it has not been created.")
        elif self.success_1 and self.success_2 and self.success_3:
            shutil.rmtree(self.tmp_dir)
