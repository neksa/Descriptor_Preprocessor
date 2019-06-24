import logging
import unittest
import os
import pickle
import shutil

from utils import extract_parser
import config

class TestIoncomExtract(unittest.TestCase):
    def setUp(self):
        ref_pname_cid_path = os.path.join(config.ROOT, 'tests', 'data',
                                          'ref', 'pname_cid_ioncom.pkl')
        with open(ref_pname_cid_path, 'rb') as file:
            self.ref_pname_cid = pickle.load(file)
        return

    def test_parse_ioncom(self):
        from utils import extract_parser
        ioncom_input_path = os.path.join(config.ROOT, 'tests', 'data',
                                          'input', 'ioncom.txt')
        self.assert_(os.path.isfile(ioncom_input_path))
        act_pname_cid_map = extract_parser.parse_ioncom(ioncom_input_path)
        self.assert_(isinstance(act_pname_cid_map, dict))
        self.assertDictEqual(self.ref_pname_cid, act_pname_cid_map)
        return

class TestPrositeExtract(unittest.TestCase):
    def setUp(self):
        ref_pname_cid_path = os.path.join(config.ROOT, 'tests', 'data', 'ref',
                                          'pname_cid_prosite.pkl')
        with open(ref_pname_cid_path, 'rb') as file:
            self.ref_pname_cid = pickle.load(file)
        return

    def test_parse_prosite(self):
        prosite_input = os.path.join(config.ROOT, 'tests', 'data', 'input',
                                          'prosite_extract.txt')
        self.assert_(os.path.isfile(prosite_input))
        act_pname_cid_map = extract_parser.parse_prosite(
            prosite_input, config.prosite_pdb_list)
        self.assert_(isinstance(act_pname_cid_map, dict))
        self.assertDictEqual(self.ref_pname_cid, act_pname_cid_map)
        return

class TestMotifFinder(unittest.TestCase):
    def setUp(self):
#         Set up test_dir
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
        self.assert_(os.path.isfile(pname_cid_ioncom_path))
        with open(pname_cid_ioncom_path, 'rb') as file:
            self.pname_cid_ioncom = pickle.load(file)

        pname_cid_prosite_path = os.path.join(
            ref_folder, 'pname_cid_prosite.pkl')
        self.assert_(os.path.isfile(pname_cid_prosite_path))
        with open(pname_cid_prosite_path, 'rb') as file:
            self.pname_cid_prosite = pickle.load(file)

        self.pdb_folder = os.path.join(config.ROOT, 'data', 'input',
                                       'pdb_files')

        ptr_props_prosite_path = os.path.join(
            ref_folder, 'ptr_props_prosite.pkl')
        with open(ptr_props_prosite_path, 'rb') as file:
            self.ref_ptr_props_prosite = pickle.load(file)

        ptr_props_ioncom_path = os.path.join(
            ref_folder, 'ptr_props_ioncom.pkl')
        with open(ptr_props_ioncom_path, 'rb') as file:
            self.ref_ptr_props_ioncom = pickle.load(file)

        self.ref_meme_path = os.path.join(ref_folder, 'meme.txt')
        return

    def test_find_motif_pos_prosite(self):
        from utils import motif_finder
        act_ptr_props = motif_finder.find_motif_pos(
            self.pname_cid_prosite, self.pdb_folder, type='prosite',
            store_dir=config.store_dir, replace_existing=True)
        self.assert_(isinstance(act_ptr_props, dict))
        self.assertDictEqual(act_ptr_props, self.ref_ptr_props_prosite)
        return True

    def test_find_motif_pos_ioncom(self):
        from utils import motif_finder
        act_ptr_props = motif_finder.find_motif_pos(
            self.pname_cid_ioncom, self.pdb_folder, type='ioncom',
            store_dir=config.store_dir, replace_existing=True,
            ref_meme_txt=self.ref_meme_path)
        self.assert_(isinstance(act_ptr_props, dict))
        self.assertDictEqual(act_ptr_props, self.ref_ptr_props_ioncom)
        return True

    def tearDown(self):
        if not os.path.isdir(self.tmp_dir):
            logging.warning(
                f"test_store_dir <{self.tmp_dir}> does not exists when "
                f"running TearDown. Check to see why it has not been created.")
        else:
            shutil.rmtree(self.tmp_dir)
        return

