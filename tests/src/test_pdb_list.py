import logging
import unittest
import os
import pickle
import shutil

from utils import extract_parser
import config

class TestDatalonExtract(unittest.TestCase):
    def setUp(self):
        ref_pname_cid_path = os.path.join(config.ROOT, 'tests', 'data',
                                          'ref', 'pname_cid_datalon.pkl')
        with open(ref_pname_cid_path, 'rb') as file:
            self.ref_pname_cid = pickle.load(file)
        return

    def test_parse_datalon(self):
        from utils import extract_parser
        datalon_input_path = os.path.join(config.ROOT, 'tests', 'data',
                                          'input', 'datalon_shortened.txt')
        self.assert_(os.path.isfile(datalon_input_path))
        act_pname_cid_map = extract_parser.parse_datalon(datalon_input_path)
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
        prosite_input_path = os.path.join(config.ROOT, 'tests', 'data', 'input',
                                          'prosite_extract_cropped.txt')
        self.assert_(os.path.isfile(prosite_input_path))
        act_pname_cid_map = extract_parser.parse_prosite_extract(
            prosite_input_path, config.prosite_pdb_list)
        self.assert_(isinstance(act_pname_cid_map, dict))
        self.assertDictEqual(self.ref_pname_cid, act_pname_cid_map)
        return

class TestPDBList(unittest.TestCase):
    def setUp(self):
#         Set up test_dir
        self.tmp_dir = os.path.join(config.ROOT, 'tests', 'data', 'tmp_PDBList')
        if os.path.isdir(self.tmp_dir):
            logging.warning(f"Test {self.tmp_dir} folder not deleted "
                            f"previously. Deleting.")
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)
        os.mkdir(os.path.join(self.tmp_dir, 'input'))

        self.store_dir = os.path.join(self.tmp_dir, 'store')
        os.mkdir(self.store_dir)

        pname_cid_datalon_path = os.path.join(config.ROOT, 'tests', 'data',
                                           'ref', 'pname_cid_datalon.pkl')
        self.assert_(os.path.isfile(pname_cid_datalon_path))
        with open(pname_cid_datalon_path, 'rb') as file:
            self.pname_cid_datalon = pickle.load(file)

        pname_cid_prosite_path = os.path.join(config.ROOT, 'tests', 'data',
                                           'ref', 'pname_cid_prosite.pkl')
        self.assert_(os.path.isfile(pname_cid_prosite_path))

        with open(pname_cid_prosite_path, 'rb') as file:
            self.pname_cid_prosite = pickle.load(file)

        # self.pdb_folder = os.path.join(self.store_dir, 'pdb_folder')
        # os.mkdir(self.pdb_folder)
        self.pdb_folder = os.path.join(config.ROOT, 'data', 'input',
                                       'pdb_files')

        ptr_props_prosite_path = os.path.join(config.ROOT, 'tests', 'data',
                                           'ref', 'ptr_props_cropped_prosite.pkl')
        with open(ptr_props_prosite_path, 'rb') as file:
            self.ref_ptr_props_prosite = pickle.load(file)

        ptr_props_datalon_path = os.path.join(config.ROOT, 'tests', 'data',
                                              'ref', 'ptr_props_cropped_datalon.pkl')
        with open(ptr_props_datalon_path, 'rb') as file:
            self.ref_ptr_props_datalon = pickle.load(file)

        self.ref_meme_path = os.path.join(config.ROOT, 'tests', 'data',
                                           'ref', 'meme.txt')
        return

    def test_get_pdb_list_prosite(self):
        from utils import pdb_list
        act_ptr_props = pdb_list.get_pdb_list(
            self.pname_cid_prosite, self.pdb_folder, type='prosite_extract',
            store_dir=config.store_dir, replace_existing=True)
        self.assert_(isinstance(act_ptr_props, dict))
        self.assertDictEqual(act_ptr_props, self.ref_ptr_props_prosite)
        return True

    def test_get_pdb_list_datalon(self):
        from utils import pdb_list
        act_ptr_props = pdb_list.get_pdb_list(
            self.pname_cid_datalon, self.pdb_folder, type='datalon',
            store_dir=config.store_dir, replace_existing=True,
            ref_meme_txt=self.ref_meme_path)
        self.assert_(isinstance(act_ptr_props, dict))
        self.assertDictEqual(act_ptr_props, self.ref_ptr_props_datalon)
        return True

    def tearDown(self):
        if not os.path.isdir(self.tmp_dir):
            logging.warning(
                f"test_store_dir <{self.tmp_dir}> does not exists when "
                f"running TearDown in TestPDBList. Check to see why it has "
                f"not been created. ")
        else:
            shutil.rmtree(self.tmp_dir)
        return

