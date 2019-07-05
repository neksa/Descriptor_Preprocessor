import logging
import os
import pickle
import shutil
import unittest

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

        self.pname_cid_ioncom_path = os.path.join(
            ref_folder, 'pname_cid_ioncom.pkl')
        self.assertTrue(os.path.isfile(self.pname_cid_ioncom_path))

        self.pname_cid_prosite_path = os.path.join(
            ref_folder, 'pname_cid_prosite.pkl')
        self.assertTrue(os.path.isfile(self.pname_cid_prosite_path))

        self.pdb_folder = os.path.join(config.ROOT, 'data', 'input',
                                       'pdb_files')

        self.motif_pos_prosite_path = os.path.join(
            ref_folder, 'motif_pos_prosite.pkl')
        self.assertTrue(os.path.isfile(self.motif_pos_prosite_path))

        self.motif_pos_ioncom_path = os.path.join(
            ref_folder, 'motif_pos_ioncom.pkl')
        self.assertTrue(os.path.isfile(self.motif_pos_ioncom_path))

        self.ref_meme_path = os.path.join(ref_folder, 'meme.txt')
        self.success = False

    def test_find_motif_pos_prosite(self):
        with open(self.pname_cid_prosite_path, 'rb') as file:
            pname_cid_prosite = pickle.load(file)
        with open(self.motif_pos_prosite_path, 'rb') as file:
            ref_motif_pos_prosite = pickle.load(file)
        act_motif_pos = motif_finder.find_motif_pos(
            pname_cid_prosite, self.pdb_folder, process='meme',
            store_dir=self.store_dir, replace_existing=True)
        self.assertIsInstance(act_motif_pos, dict)
        self.assertDictEqual(act_motif_pos, ref_motif_pos_prosite)
        self.success = True

    def test_find_motif_pos_ioncom(self):
        with open(self.pname_cid_ioncom_path, 'rb') as file:
            pname_cid_ioncom = pickle.load(file)
        with open(self.motif_pos_ioncom_path, 'rb') as file:
            ref_motif_pos_ioncom = pickle.load(file)
        act_motif_pos = motif_finder.find_motif_pos(
            pname_cid_ioncom, self.pdb_folder, process='mast',
            store_dir=self.store_dir, replace_existing=True,
            ref_meme_txt=self.ref_meme_path)
        self.assertIsInstance(act_motif_pos, dict)
        self.assertDictEqual(act_motif_pos, ref_motif_pos_ioncom)
        self.success = True

    def test_find_motif_pos_prosite_no_pdb_available(self):
        with open(self.pname_cid_prosite_path, 'rb') as file:
            pname_cid_prosite = pickle.load(file)

        empty_pdb_folder = os.path.join(self.tmp_dir, 'empty_pdb_files')
        if os.path.isdir(empty_pdb_folder):
            logging.warning(f"Empty_pdb_folder in <{empty_pdb_folder}> not "
                            f"deleted. Deleting.")
            shutil.rmtree(empty_pdb_folder)
        os.mkdir(empty_pdb_folder)
        shortened_pname_cid = dict()

        for i, key in enumerate(pname_cid_prosite.keys()):
            if i == 2:
                break
            shortened_pname_cid[key] = pname_cid_prosite[key]
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
        #     todo: check if .pdb file is legit?
        self.success = True

    def tearDown(self):
        if not os.path.isdir(self.tmp_dir):
            logging.warning(
                f"test_store_dir <{self.tmp_dir}> does not exists when "
                f"running TearDown. Check to see why it has not been created.")
        elif self.success:
            shutil.rmtree(self.tmp_dir)
